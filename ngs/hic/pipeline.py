"""Pipeline to deal with Hi-C data"""

# import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import rpy2.robjects as robjects
import subprocess


def run_hicpro(
        input_path, output_path,
        hicpro_path='/home/wangyn/software/HiC-Pro_2.7.8/bin/HiC-Pro',
        config_template_path='/home/wangyn/lustre/pipeline/Hi-C/config_template_mboi.txt'):
    """Run HicPro on knight server

    Args:
        input_path: string, the path of input files
        output_path: string, the path of output
        hicpro_path: string, install path of HiC-Pro,
        has been set for knight server
        config_template_path: string, the path of config template,
        has been set for knight server
    Returns:
        0
    """
    print('Precessing data in {}, the outputs are in {}'.format(input_path,
                                                                output_path))
    cmd = 'nohup {} -i {} -o {} -c {} &'.format(hicpro_path, input_path,
                                                output_path,
                                                config_template_path)
    print('Run: {}'.format(cmd))
    subprocess.Popen(cmd, shell=True)
    return 0


def copy_to_db(df, table_name,
               engine=sqlalchemy.create_engine(
                   'postgresql://wangyn@localhost/lung_cancer_db'),
               if_exists='fail'):
    """Copy a DataFrame to database

    Args:
        df: the DataFrame to be copied
        table_name: the name to be assigned in the database
        engine: sqlalchemy engine
        if_exists: what to do if the table has existed
    Returns:
        0
    """
    try:
        df.to_sql(table_name, engine, if_exists=if_exists, index=False)
        return 0
    except ValueError:
        print('The table {} has existed.'.format(table_name))
        return 1


def load_hicpro_outputs(matrix_path, bed_path):
    """Load the outputs of HiC-Pro

    Args:
        matrix_path: the path of matrix file
        bed_path: the path of bed file
    Returns:
        hicpro_dict, a dict containing both matrix and bed data outputed from HiC-Pro
    """
    hicpro_dict = {}
    matrix_data = pd.read_table(matrix_path, header=None,
                                names=['bin_1', 'bin_2', 'counts'])
    bed_data = pd.read_table(bed_path, header=None,
                             names=['chrom', 'bin_start', 'bin_end',
                                    'bin_index'])
    hicpro_dict['matrix'] = matrix_data
    hicpro_dict['bed'] = bed_data
    return hicpro_dict


def load_hic(matrix_path, bed_path):
    """Load the outputs of HiC-Pro

    Args:
        matrix_path: the path of matrix file
        bed_path: the path of bed file
    Returns:
        The full matrix of Hi-C data
    """
    matrix_data = pd.read_table(matrix_path, header=None,
                                names=['bin_1', 'bin_2', 'counts'])
    bed_data = pd.read_table(bed_path, header=None,
                             names=['chrom', 'bin_start', 'bin_end',
                                    'bin_index'])
    chroms = list(set(bed_data.chrom))
    hic_full_matrix = {}
    for chrom in chroms:
        print('Processing {}...'.format(chrom))
        sub_bed = bed_data[bed_data.chrom == chrom]
        sub_matrix = matrix_data[matrix_data.bin_1.isin(sub_bed.bin_index) &
                                 matrix_data.bin_2.isin(sub_bed.bin_index)]
        sub_full_mat = pd.DataFrame(
            np.zeros([sub_bed.shape[0], sub_bed.shape[0] + 3]),
            index=list(sub_bed.bin_index),
            columns=['chromosome', 'bin_start', 'bin_end'] + list(
                sub_bed.bin_index))
        sub_full_mat.iloc[:, 0] = chrom
        sub_full_mat.iloc[:, 1] = list(sub_bed.bin_start)
        sub_full_mat.iloc[:, 2] = list(sub_bed.bin_end)
        bin_1 = list(sub_matrix.bin_1)
        bin_2 = list(sub_matrix.bin_2)
        for i in range(sub_matrix.shape[0]):
            index_1 = bin_1[i]
            index_2 = bin_2[i]
            sub_full_mat.loc[index_1, index_2] += 1
            if index_1 != index_2:
                sub_full_mat.loc[index_2, index_1] += 1
        hic_full_matrix[chrom] = sub_full_mat
    print('Finish.')
    return hic_full_matrix


def call_tad(matrix_file, output_file, window_size=5):
    """Call TAD by TopDom

    Args:
        matrix_file: the path of matrix file
        output_file: the path of output, the outputs are .bed, .domain
        and .binSignal
        window_size: window size, For more details, see TopDom's manual:
                     http://zhoulab.usc.edu/TopDom/TopDom%20Manual_v0.0.2.pdf
    Returns:
        0
    """

    rTopDom = robjects.r("""
                         source('/home/wangyn/lustre/pipeline/Hi-C/src/TopDom_v0.0.2.R')
                         TopDom
                         """)

    matrix_file = robjects.StrVector([matrix_file])
    window_size = robjects.IntVector([window_size])
    out_file = robjects.StrVector([output_file])
    rTopDom(matrix_file, window_size, out_file)
    return 0


def hic_pipeline(sample_name, input_dir, output_dir, resolution, window_size=5):
    """This function use Hic-Pro's outputs as inputs, automatically generate
       full matrix from sparse matrix and call TAD for each chromosome.
       For each sample or Hi-C experiment, the function will generate a project
       which is a directory containing the files for full matrix and TADs.

       Args:
           sample_name: the name of the experiment
           input_dir: results of Hic-Pro
           output_dir: the output path
           resolution: the resolution of data
           window_size: the window size used when call TAD, default is 5
        Returns:
            0
    """

    print('Run Hi-C pipeline...')
    output_path = os.path.join(output_dir, sample_name)
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    else:
        print('The dir {} exists. Please choose another path.'.format(
            output_path))
        return 1

    # generate full matrix
    hic_matrix_path = os.path.join(output_path, 'hic_matrix')
    raw_out = os.path.join(hic_matrix_path, 'raw')
    iced_out = os.path.join(hic_matrix_path, 'iced')
    os.mkdir(hic_matrix_path)
    os.mkdir(raw_out)
    os.mkdir(iced_out)

    # raw
    print('Generate raw Hi-C matrix...')
    raw_dir = os.path.join(input_dir, 'raw', str(resolution))
    raw_files = os.listdir(raw_dir)
    raw_files.sort()
    raw_matrix_path = os.path.join(raw_dir, raw_files[0])
    bed_path = os.path.join(raw_dir, raw_files[1])

    raw_mat_dict = load_hic(raw_matrix_path, bed_path)
    for key in raw_mat_dict.keys():
        tmp = raw_mat_dict[key]
        tmp.to_csv(os.path.join(raw_out, key), sep='\t', index=False,
                   header=False)

    # ice
    print('Generate iced Hi-C matrix...')
    iced_dir = os.path.join(input_dir, 'iced', str(resolution))
    iced_files = os.listdir(iced_dir)
    iced_matrix_path = os.path.join(iced_dir, iced_files[0])

    iced_mat_dict = load_hic(iced_matrix_path, bed_path)
    for key in iced_mat_dict.keys():
        tmp = iced_mat_dict[key]
        tmp.to_csv(os.path.join(iced_out, key), sep='\t', index=False,
                   header=False)

    # Call TAD
    print('Call TAD...')
    fail_files = []
    tad_out = os.path.join(output_path, 'tad')
    os.mkdir(tad_out)
    for key in iced_mat_dict.keys():
        matrix_file = os.path.join(iced_out, key)
        output_file = os.path.join(tad_out, key)
        try:
            call_tad(matrix_file, output_file, window_size=window_size)
        except:
            fail_files.append(matrix_file)
    print('Fail to call TAD from file: {}'.format(matrix_file))
    return 0
