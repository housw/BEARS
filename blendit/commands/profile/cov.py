# -*- coding: utf-8 -*-

import os
import logging
import numpy as np
import pandas as pd
import multiprocessing
from blendit.utils.external import run_bamcov
from blendit.utils.common import normalizer
from blendit.utils.common import CommandException
from blendit.utils.common import folder_exists
from blendit.utils.common import create_directory
from blendit.utils.common import emit_file_exist_warning


_logger = logging.getLogger("BlendIt")


def bamcov_worker(param_dict):
    pid = os.getpid()
    try:
        bamcov_file = os.path.join(param_dict['output_dir'], param_dict['depth_prefix'] + ".bamcov")
        try:
            bamcov_file = emit_file_exist_warning(filename=bamcov_file, force=param_dict['force'])
        except Exception as e:
            _logger.info(e)
            bamcov_file = run_bamcov(**param_dict)
            _logger.info("bamcov has finished counting {0} on pid {1}".format(param_dict['bam_file'], pid))
        return bamcov_file
    except Exception as e:
        err_msg = "bamcov on pid {0} failed due to {1}".format(pid, e)
        _logger.error(err_msg)
        raise CommandException(err_msg)


def parallel_bamcov(bam_file_list, output_dir, min_read_len=30, min_MQ=0, min_BQ=0, cpus=10, force=False):

    # use dict as kwargs of run_bamcov
    param_dicts = []
    sample_names = []
    for bam_file in bam_file_list:
        if os.path.isfile(bam_file):
            basename = os.path.basename(bam_file)
            depth_prefix = os.path.splitext(basename)[0]
            param_dict = {'bam_file': bam_file,
                          'depth_prefix': depth_prefix,
                          'output_dir': output_dir,
                          'min_read_len': min_read_len,
                          'min_MQ': min_MQ,
                          'min_BQ': min_BQ,
                          'force': force}
            param_dicts.append(param_dict)
            sample_names.append(depth_prefix)
        else:
            err_msg = "{0} is not a file, please make sure you have correct input!".format(bam_file)
            _logger.error(err_msg)
            raise CommandException(err_msg)

    pool = multiprocessing.Pool(cpus)
    depth_files = pool.map(bamcov_worker, param_dicts)
    pool.close()

    return sample_names, depth_files


def write_length_and_depth_file(sample_names, depth_files, output_length_file, output_depth_file):

    # write contig length file
    first_df = pd.read_csv(depth_files[0], sep='\t', header=0, index_col=0)
    length_df = first_df[['endpos']].copy()
    length_df.rename(columns={'endpos':'Length'}, inplace=True)
    length_df.index.name = "Contig_ID"
    length_df.to_csv(output_length_file, sep='\t', header=True, index=True)

    # merge depth profiles and write depth file
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

    return length_df, depth_df


def parallel_calculate_contig_depth_from_bam_files(input_bam_file_list, output_dir, prefix, min_read_len=30,
                                                   min_MQ=0, min_BQ=0, cpus=10, force=False,
                                                   scale_func=np.log10, read_length=250):
    # output file handling
    if not folder_exists(output_dir):
        create_directory(output_dir)
    output_length_file = os.path.join(output_dir, prefix + "_contig_length.tsv")
    output_depth_file = os.path.join(output_dir, prefix + "_coverage.tsv")
    output_depth_prior = os.path.join(output_dir, prefix + "_coverage_prior.tsv")
    output_depth_norm = os.path.join(output_dir, prefix + "_coverage_prior_norm.tsv")

    _logger.info("parallel calculating contig coverage using bamcov ...")
    try:
        emit_file_exist_warning(filename=output_depth_file, force=force)
    except Exception as e:
        _logger.info(e)
        sample_names, depth_files = parallel_bamcov(input_bam_file_list, output_dir, min_read_len=min_read_len,
                                                    min_MQ=min_MQ, min_BQ=min_BQ, cpus=cpus, force=force)
        length_df, depth_df = write_length_and_depth_file(sample_names, depth_files, output_length_file, output_depth_file)

    length_df.sort_index(inplace=True)
    depth_df.sort_index(inplace=True)

    # add one read to each contig as prior
    _logger.info("spiking one aligned read to each contig ...")
    try:
        emit_file_exist_warning(filename=output_depth_prior, force=force)
    except Exception as e:
        _logger.info(e)
        prior_df = read_length / length_df.Length
        for col in depth_df.columns:
            depth_df[col] += prior_df
        depth_df.to_csv(output_depth_prior, sep="\t", header=True, index=True)

    # log-scale transform the data, run normalization
    _logger.info("normalizing depth profiles ...")
    try:
        emit_file_exist_warning(filename=output_depth_norm, force=force)
    except Exception as e:
        _logger.info(e)
        output_depth_norm = normalizer(input_freq_file=output_depth_prior, output_file=output_depth_norm,
                                       scale_func=scale_func)

    return output_length_file, output_depth_norm


def calculate_contig_depth_from_bam_files(input_bam_file_list, output_dir, prefix,
                                          min_read_len=30, min_MQ=0, min_BQ=0, force=False,
                                          scale_func=np.log10, read_length=250):
    # output file handling
    if not folder_exists(output_dir):
        create_directory(output_dir)
    output_length_file = os.path.join(output_dir, prefix + "_contig_length.tsv")
    output_depth_file = os.path.join(output_dir, prefix + "_coverage.tsv")
    output_depth_prior = os.path.join(output_dir, prefix + "_coverage_prior.tsv")
    output_depth_norm = os.path.join(output_dir, prefix + "_coverage_prior_norm.tsv")

    # run bamcov iteratively
    depth_files = []
    sample_names = []
    for bam_file in input_bam_file_list:
        if os.path.isfile(bam_file):
            basename = os.path.basename(bam_file)
            depth_prefix = os.path.splitext(basename)[0]
            bamcov_file = os.path.join(output_dir, depth_prefix + ".bamcov")
            try:
                emit_file_exist_warning(filename=bamcov_file, force=force)
            except Exception as e:
                _logger.info(e)
                bamcov_file = run_bamcov(bam_file, depth_prefix, output_dir,
                                         min_read_len=min_read_len, min_MQ=min_MQ, min_BQ=min_BQ, force=force)
            sample_names.append(depth_prefix)
            depth_files.append(bamcov_file)
        else:
            err_msg = "{0} is not a bam file, please make sure you have correct input!".format(bam_file)
            _logger.error(err_msg)
            raise CommandException(err_msg)

    # write contig length file
    _logger.info("generating contig length file ...")
    first_df = pd.read_csv(depth_files[0], sep='\t', header=0, index_col=0)
    length_df = first_df[['endpos']].copy()
    length_df.rename(columns={'endpos':'Length'}, inplace=True)
    length_df.index.name = "Contig_ID"
    length_df.to_csv(output_length_file, sep='\t', header=True, index=True)

    # merge depth profiles and write depth file
    _logger.info("concatenating coverage from all the samples ...")
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
    depth_df = pd.concat(all_depth_dfs, axis=1, sort=True)
    depth_df.index.name="Contig_ID"
    depth_df.to_csv(output_depth_file, sep="\t", header=True, index=True)

    length_df.sort_index(inplace=True)
    depth_df.sort_index(inplace=True)

    # add one read to each contig as prior
    prior_df = read_length / length_df.Length
    for col in depth_df.columns:
        depth_df[col] += prior_df
    depth_df.to_csv(output_depth_prior, sep="\t", header=True, index=True)

    # log-scale transform the data, run normalization
    norm_depth_file = normalizer(input_freq_file=output_depth_prior, output_file=output_depth_norm, scale_func=scale_func)

    return output_length_file, output_depth_norm
