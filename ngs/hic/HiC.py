"""HiC class, integrate the functions in pipeline into one class"""
import matplotlib.gridspec as grd
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import rpy2.robjects as robjects
import plotly.plotly as py
import plotly.graph_objs as go


class HiC:
    """Results of HiC experiment"""

    def __init__(self, chrom, resolution, full_matrix=None,
                 matrix=None, bed=None, tad=None, compartment=None,
                 insulation_score=None, directionality_index=None):
        """Initialization

        Args:
            full_matrix: Hi-C matrix,
            chrom: the chromosome that described by the matrix
            resolution: resolution of Hi-C matrix

        Returns:
            None
        """
        self.full_matrix = full_matrix
        self.chrom = chrom
        self.resolution = resolution
        self.matrix = matrix
        self.bed = bed
        self.tad = tad
        self.compartment = compartment
        self.insulation_score = insulation_score
        self.directionality_index = directionality_index

    def load_hicpro_outputs(self, matrix_path, bed_path):
        """Load HiC-Pro outputs (matrix and bed) from files

        Args:
            matrix_path: path of .matrix file
            bed_path: path of .bed file

        Returns:
            None
        """
        self.matrix = pd.read_table(matrix_path, header=None,
                                    names=['bin_1', 'bin_2', 'counts'])
        self.bed = pd.read_table(bed_path, header=None,
                                 names=['chrom', 'bin_start', 'bin_end',
                                        'bin_index'])

    def load_from_sql(self, matrix_name, bed_name, engine, sub_set=True):
        """Load HiC-Pro outputs (matrix and bed) from sql

        Args:
            matrix_name: table name of .matrix file
            bed_name: table name of .bed file
            engine: sqlalchemy engine
            sub_set: only retain the data corresponding to self.chrom

        Returns:
            None
        """
        self.matrix = pd.read_sql(matrix_name, engine)
        self.bed = pd.read_sql(bed_name, engine)
        if sub_set:
            assert self.chrom is not None
            self.bed = self.bed[self.bed.chrom == self.chrom]
            self.matrix = self.matrix[
                self.matrix.bin_1.isin(self.bed.bin_index) &
                self.matrix.bin_2.isin(self.bed.bin_index)]

    def get_full_mat(self):
        """Transform sparse matrix to full matrix

        Returns:
            None
        """
        full_mat = pd.DataFrame(
            np.zeros([self.bed.shape[0], self.bed.shape[0] + 3]),
            index=list(self.bed.bin_index),
            columns=['chromosome', 'bin_start', 'bin_end'] + list(
                self.bed.bin_index))
        full_mat.iloc[:, 0] = self.chrom
        full_mat.iloc[:, 1] = list(self.bed.bin_start)
        full_mat.iloc[:, 2] = list(self.bed.bin_end)
        bin_1 = list(self.matrix.bin_1)
        bin_2 = list(self.matrix.bin_2)
        for i in range(self.matrix.shape[0]):
            if i % 5000 == 0:
                print('---------{}%'.format(
                    round(i / self.matrix.shape[0] * 100, 2)))
            index_1 = bin_1[i]
            index_2 = bin_2[i]
            full_mat.loc[index_1, index_2] += self.matrix.iloc[i, 2]
            if index_1 != index_2:
                full_mat.loc[index_2, index_1] += self.matrix.iloc[i, 2]
        self.full_matrix = full_mat

    def get_tad(self, window_size=5, matrix_file='./tmp_matrix_file',
                output='./tmp_output'):
        """Call TADs from Hi-C matrix, the algorithm used here is TopDom,
        For more details, see TopDom's manual:
        http://zhoulab.usc.edu/TopDom/TopDom%20Manual_v0.0.2.pdf

        Args:
            window_size: the window size of TopDom
            matrix_file: the path of matrix file
            output: the path of output
        Returns:
            None
        """
        self.full_matrix.to_csv(matrix_file, sep='\t', index=False,
                                header=False)
        rTopDom = robjects.r("""
        source('/home/wangyn/lustre/pipeline/Hi-C/src/TopDom_v0.0.2.R')
        TopDom
        """)
        r_matrix_file = robjects.StrVector([matrix_file])
        r_window_size = robjects.IntVector([window_size])
        r_out_file = robjects.StrVector([output])
        rTopDom(r_matrix_file, r_window_size, r_out_file)
        tad = pd.read_table('./tmp_output.domain')
        os.remove(matrix_file)
        os.remove(output + '.domain')
        os.remove(output + '.bed')
        os.remove(output + '.binSignal')
        self.tad = tad

    def get_compartment(self):
        pass

    def get_insulation_score(self, window_size=500e3):
        """Compute insulation score.
        For more deatails about insulation score, see:
        http://www.ncbi.nlm.nih.gov/pubmed/26030525

        Args:
            window_size: the size of window sliding along matrix's diagonal

        Returns:
            None
        """
        assert self.full_matrix is not None, 'The matrix is None'
        bin_size = window_size / self.resolution
        assert bin_size % 1 == 0, 'Bin size is not an integer'
        bin_size = int(bin_size)
        chrom_len = self.full_matrix.shape[0]
        insulation_score = pd.DataFrame(np.zeros([chrom_len - 2 * bin_size, 3]))
        insulation_score.columns = ['bin_start', 'bin_end', 'insulation_score']
        for i in range(bin_size, chrom_len - bin_size):
            iscore = np.sum(np.sum(self.full_matrix.iloc[(i - bin_size):i,
                                   (i + 4):(i + bin_size + 4)]))
            insulation_score.iloc[i - bin_size, 0] = self.full_matrix.iloc[i, 1]
            insulation_score.iloc[i - bin_size, 1] = self.full_matrix.iloc[i, 2]
            insulation_score.iloc[i - bin_size, 2] = iscore
        iscore_all = insulation_score.loc[:, 'insulation_score']
        insulation_score.loc[:, 'insulation_score'] = np.log2(
            iscore_all / np.mean(iscore_all))
        self.insulation_score = insulation_score

    def get_directionality_index(self, upstream=2e6, downstream=2e6):
        """Compute directionality index.
        For more details about directionality index, see:
        http://www.ncbi.nlm.nih.gov/pubmed/22495300

        Args:
            upstream: the upstream of bps to compute A
            downstream: the downstream of bps to compute B

        Returns:
            None
        """
        assert self.full_matrix is not None, 'The matrix is None'
        chrom_len = self.full_matrix.shape[0]
        upstream_size = upstream / self.resolution
        downstream_size = downstream / self.resolution
        assert upstream_size % 1 == 0, 'Fail to exact division'
        assert downstream_size % 1 == 0, 'Fail to exact division'
        upstream_size = int(upstream_size)
        downstream_size = int(downstream_size)
        di = np.zeros(chrom_len)
        for i in range(chrom_len):
            up_start = np.max([0, i - upstream_size])
            down_end = np.min([chrom_len, i + downstream_size])
            a = np.sum(self.full_matrix.iloc[i, (up_start + 3):(i + 3)])
            # don't include the bin itself
            b = np.sum(self.full_matrix.iloc[i, (i + 4):(down_end + 4)])
            e = (a + b) / 2
            di[i] = ((b - a) / np.abs(b - a)) * \
                    ((np.square((a - e)) / e) + (np.square((b - e)) / e))
        directionality_index = pd.DataFrame(np.zeros([chrom_len, 3]))
        directionality_index.columns = ['bin_start', 'bin_end', 'di']
        directionality_index.loc[:, 'bin_start'] = list(
            self.full_matrix.loc[:, 'bin_start'])
        directionality_index.loc[:, 'bin_end'] = list(
            self.full_matrix.loc[:, 'bin_end'])
        directionality_index.loc[:, 'di'] = di
        self.directionality_index = directionality_index

    def plot(self, max_percentile=95, add_tad=False, add_di=False,
             add_is=False):
        """Visualize Hi-C data by pyplot

        Args:
            max_percentile: the contact number above this percentile will be
                reduced to this number
            add_tad: if add TAD to the plot
            add_di: if add directionality index to the plot
            add_is: if add insulation score to the plot

        Returns:
            None
        """
        subplot_num = 1
        height_ratios = [12]
        width_ratios = [10, 1]
        if add_di:
            subplot_num += 1
            height_ratios.append(3)
        if add_is:
            subplot_num += 1
            height_ratios.append(3)

        # log2 Hi-C matrix and compress big number
        log2_mat = np.log2(self.full_matrix.iloc[:, 3:].values + 1)
        max_val = np.percentile(log2_mat, max_percentile)
        log2_mat[log2_mat > max_val] = max_val

        # initialization
        fig = plt.figure(figsize=[15, 15])
        gs = grd.GridSpec(subplot_num, 2, height_ratios=height_ratios,
                          width_ratios=width_ratios, wspace=0.1)

        # heat map
        ax1 = plt.subplot(gs[0])
        mesh = ax1.pcolorfast(log2_mat, cmap=plt.cm.Reds, norm=None)
        ax1.get_xaxis().set_ticks([])
        ax1.get_yaxis().set_ticks([])
        ax1.invert_yaxis()

        # tad
        if add_tad:
            assert self.tad is not None, 'Please run get_tad first.'
            for i in range(self.tad.shape[0]):
                if self.tad['tag'][i] == 'domain':
                    from_id = self.tad['from.id'][i]
                    to_id = self.tad['to.id'][i]
                    ax1.hlines(from_id, xmin=from_id, xmax=to_id,
                               linestyle='--', linewidth=3, colors='blue')
                    ax1.vlines(to_id, ymin=from_id, ymax=to_id,
                               linestyle='--', linewidth=3, colors='blue')
                    ax1.vlines(from_id, ymin=from_id, ymax=to_id,
                               linestyle='--', linewidth=3, colors='blue')
                    ax1.hlines(to_id, xmin=from_id, xmax=to_id,
                               linestyle='--', linewidth=3, colors='blue')

        # color bar for heat map
        colorAx = plt.subplot(gs[1])
        ax2 = plt.colorbar(mesh, cax=colorAx)

        # directionality index
        if add_di:
            di_x_vals = self.directionality_index.bin_start.values
            di_y_vals = self.directionality_index.di.values
            ax3 = plt.subplot(gs[2])
            ax3.fill_between(di_x_vals, 0, di_y_vals, where=di_y_vals >= 0,
                             facecolor='orange', interpolate=True)
            ax3.fill_between(di_x_vals, 0, di_y_vals, where=di_y_vals < 0,
                             facecolor='skyblue', interpolate=True)
            ax3.set_xlim([np.min(self.full_matrix.bin_start.values),
                          np.max(self.full_matrix.bin_start.values)])
            ax3.set_ylabel('DI')

        # insulation score
        if add_is:
            is_x_vals = self.insulation_score.bin_start.values
            is_y_vals = self.insulation_score.insulation_score.values
            ax4 = plt.subplot(gs[4])
            ax4.fill_between(is_x_vals, 0, is_y_vals, where=is_y_vals >= 0,
                             facecolor='palegoldenrod', interpolate=True)
            ax4.fill_between(is_x_vals, 0, is_y_vals, where=is_y_vals < 0,
                             facecolor='plum', interpolate=True)
            ax4.set_xlim([np.min(self.full_matrix.bin_start.values),
                          np.max(self.full_matrix.bin_start.values)])
            ax4.set_ylabel('IS')

        plt.show()

    def iplot(self, max_percentile=95):
        """Plot Hi-C data by plotly

        Args:
            max_percentile: the contact number above this percentile will be
                reduced to this number

        Returns:
            None
        """
        log2_mat = np.log2(self.full_matrix.iloc[:, 3:].values + 1)
        max_val = np.percentile(log2_mat, max_percentile)
        log2_mat[log2_mat > max_val] = max_val
        log2_mat = log2_mat.tolist()
        hm_data = [go.Heatmap(z=log2_mat, colorscale=['White', 'Red'])]
        fig = go.Figure(data=hm_data)
        fig['layout'].update(
            width=700,
            height=700,
            autosize=False
        )
        return py.iplot(fig)
