"""Pipelines to deal with ChIP-seq data,
mainly bed files downloaded from ENCODE"""

import numpy as np
import pandas as pd
import sqlalchemy


def load_data(bed_path):
    """Load data from a bed file downloaded from ENCODE

    Args:
        bed_path: file's path
    Returns:
        df, pandas DataFrame of the bed file
    """
    return pd.read_table(bed_path, header=None)


def compute_rpm(df, resolution):
    """Compute RPM for each genomic bin

    Args:
        df: the DataFrame storing data from ENCODE bed files, outputs of load_data
        resolution: genomic bin size
    Returns:
        chipseq_dict, keys are chromosomes, values are RPM and bin_start
    """
    chroms = list(set(df.iloc[:, 0]))
    print('Detect chromosomes: {}'.format(chroms))

    chipseq_dict = {}
    for chrom in chroms:
        print('Processing {}...'.format(chrom))
        df_chrom = df[df.iloc[:, 0] == chrom]
        bin_start = (
        np.floor(df_chrom.iloc[:, 1].values / resolution) * resolution).astype(
            np.int64)
        df_chrom.loc[:, 'bin_start'] = bin_start
        df_chrom_counts = df_chrom.groupby(
            ['bin_start']).count().reset_index().iloc[:, :2]
        df_chrom_counts.columns = ['bin_start', 'rpm']
        df_chrom_counts.loc[:, 'rpm'] = df_chrom_counts.loc[:, 'rpm'] / sum(
            df_chrom_counts.loc[:, 'rpm']) * 1e6
        chipseq_dict[chrom] = df_chrom_counts
    return chipseq_dict


def save_to_sql(chipseq_dict, prefix, engine, if_exists='fail'):
    """Save the results in PosegreSQL

    Args:
        chipseq_dict: the output of test_compute_rpm, a dict store ChIP-seq RPM data
        prefix: the prefix of table name
        engine: sqlalchemy engine
        if_exists: deal with the situation that the table exists,
        default is 'fail'
    Returns:
        0
    """
    for key in chipseq_dict.keys():
        table_name = prefix + '_' + key
        print('Save data into {}'.format(table_name))
        try:
            chipseq_dict[key].to_sql(table_name, engine, if_exists=if_exists)
        except ValueError:
            print('WARNING: Table {} already exists.'.format(table_name))
    return 0


def pipeline_chipseq(bed_path, resolution, prefix,
                     engine=sqlalchemy.create_engine(
                         'postgresql://wangyn@localhost/markov3d_db'),
                     if_exists='fail'):
    """Load ChIP-seq data, compute RPM and store into database

    Args:
       bed_path: file's path
       resolution: genomic bin size
       prefix: the prefix of table name
       engine: sqlalchemy engine
       if_exists: deal with the situation that the table exists,
       default is 'fail'
    Returns:
        0
    """
    # load data
    chipseq = load_data(bed_path)

    # compute RPM
    chipseq_dict = compute_rpm(chipseq, resolution)

    # copy to database
    save_to_sql(chipseq_dict, prefix, engine, if_exists=if_exists)

    return 0
