# -*- coding: utf-8 -*-

import os
import logging
import numpy as np
import pandas as pd
from MulticoreTSNE import MulticoreTSNE as TSNE
from blendit.utils.external import run_bamcov
from blendit.utils.common import normalizer


_logger = logging.getLogger("BlendIt")


def calculate_contig_depth_from_bam_files(input_bam_file_list, output_dir, output_prefix,
                                          min_read_len=30, min_MQ=0, min_BQ=0, force=False):

    depth_files = []
    sample_names = []

    # run bamcov iteratively
    for bam_file in input_bam_file_list:
        basename = os.path.basename(bam_file)
        depth_prefix = os.path.splitext(basename)[0]
        bamcov_file = run_bamcov(bam_file, depth_prefix, output_dir,
                                 min_read_len=min_read_len, min_MQ=min_MQ, min_BQ=min_BQ, force=force)
        sample_names.append(depth_prefix)
        depth_files.append(bamcov_file)

    # write contig length file
    output_length_file = os.path.join(output_dir, output_prefix+"_contig_length.tsv")
    first_df = pd.read_csv(depth_files[0], sep='\t', header=0, index_col=0)
    length_df = first_df[['endpos']].copy()
    length_df.rename(columns={'endpos':'Length'}, inplace=True)
    length_df.index.name = "Contig_ID"
    length_df.to_csv(output_length_file, sep='\t', header=True, index=True)

    # merge depth profiles and write depth file
    output_depth_file = os.path.join(output_dir, output_prefix+".tsv")
    all_depth_dfs = []
    shape = None
    for i, depth_file in enumerate(depth_files):
        curr_depth_df = pd.read_csv(depth_file, sep='\t', header=0, index_col=0).meandepth.rename(sample_names[i])
        if not shape:
            shape = curr_depth_df.shape
            all_depth_dfs.append(curr_depth_df)
        else:
            if curr_depth_df.shape == shape:
                all_depth_dfs.append(curr_depth_df)
            else:
                sample_id = sample_names[i]
                _logger.error("sample {0} has a different dimension, this sample is temporarily ignored, "
                              "please re-generate the bam and bamcov files if you want to include it".format(sample_id))
    _logger.info("concatenating coverage from all the samples ...")
    depth_df = pd.concat(all_depth_dfs, axis=1, sort=True)
    depth_df.index.name="Contig_ID"
    depth_df.to_csv(output_depth_file, sep="\t", header=True, index=True)

    return output_length_file, output_depth_file


def normalize_contig_depth(input_length_file, input_depth_file, output_dir, output_prefix, scale_func=np.log10, read_length=250):
    """ 1) add one read to each contig as prior
        2) scale the data with numpy function
        3) normalize the data using Normalizer
    """
    length_df = pd.read_csv(input_length_file, sep="\t", index_col=0, header=0)
    length_df.sort_index(inplace=True)
    depth_df = pd.read_csv(input_depth_file, sep="\t", index_col=0, header=0)
    depth_df.sort_index(inplace=True)

    # add 1 read prior to depth_df,
    prior_depth_file  = os.path.join(output_dir, output_prefix + ".prior")
    prior_df = read_length / length_df.Length
    for col in depth_df.columns:
        depth_df[col] += prior_df
    depth_df.to_csv(prior_depth_file, sep="\t", header=True, index=True)

    # log-scale transform the data, run normalization
    norm_depth_file = normalizer(prior_depth_file, output_dir, output_prefix, scale_func=scale_func)

    return norm_depth_file
