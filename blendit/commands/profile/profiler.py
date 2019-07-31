# -*- coding: utf-8 -*-

import os
import logging
import warnings
warnings.filterwarnings(action='ignore', category=DeprecationWarning)
import pandas as pd
from .kmer import get_kmer_frequencies_for_contigs
from .codon import get_codon_frequencies_for_contigs
from .cov import parallel_calculate_contig_depth_from_bam_files
from blendit.utils.embedding import compute_embeddings
from blendit.utils.common import folder_exists
from blendit.utils.common import create_directory
from blendit.utils.common import emit_file_exist_warning


_logger = logging.getLogger("BlendIt")


def blendit_profiler(assembly, bam_files, prefix, output_dir, kmer_size, kmerfreq_scale_func,
                     codonfreq_scale_func, genetic_code, cov_scale_func, min_read_len, min_MQ, min_BQ, threads,
                     dimred, dimensions, components, read_length, force,
                     *args, **kwargs):

    # get feature profiles
    norm_kmer_file, norm_codon_file, contig_length_file, depth_file = blendit_get_profiles(assembly, bam_files,
                     prefix, output_dir, kmer_size, kmerfreq_scale_func, codonfreq_scale_func, genetic_code,
                     cov_scale_func, min_read_len, min_MQ, min_BQ, threads, read_length, force)

    # get merged feature profiles
    feature_files = blendit_merge_profiles(norm_kmer_file, norm_codon_file, depth_file, output_dir, prefix, force, threads,
                           dimred, dimensions, components)

    return feature_files


def blendit_merge_profiles(kmerfreq_file, codonfreq_file, depth_file, output_dir, prefix, force, threads,
                           dimred, dimensions, components):

    _logger.info("combining kmerfreq, codonfreq and depth profiles ...")
    merged_profile = combine_feature_tables([kmerfreq_file, codonfreq_file, depth_file],
                                            output_dir=output_dir, prefix=prefix, force=force)

    _logger.info("computing embeddings based on merged profiles...")
    feature_files = compute_embeddings(merged_profile=merged_profile, output_dir=output_dir, prefix=prefix,
                    threads=threads, n_components=dimensions, pca_components=components, dimred=dimred, force=force)

    return feature_files


def blendit_get_profiles(assembly, bam_files, prefix, output_dir, kmer_size, kmerfreq_scale_func,
                     codonfreq_scale_func, genetic_code, cov_scale_func, min_read_len, min_MQ, min_BQ, threads,
                     read_length, force, *args, **kwargs):

    # run kmer frequency calculation
    norm_kmer_file = get_kmer_frequencies_for_contigs(input_contig_file=assembly,
                                                      output_dir=os.path.join(output_dir, "kmer"),
                                                      prefix=prefix, k=kmer_size, cpus=threads,
                                                      force=force, scale_func=kmerfreq_scale_func)
    # run codon frequency calculation
    norm_codon_file = get_codon_frequencies_for_contigs(input_contig_file=assembly,
                                                        output_dir=os.path.join(output_dir, "codon"),
                                                        prefix=prefix, genetic_code=genetic_code, cpus=threads,
                                                        force=force, scale_func=codonfreq_scale_func)

    # run coverage calculation
    contig_length_file, depth_file = parallel_calculate_contig_depth_from_bam_files(input_bam_file_list=bam_files,
                                                        output_dir=os.path.join(output_dir, "cov"),
                                                        prefix=prefix, min_read_len=min_read_len,
                                                       min_MQ=min_MQ, min_BQ=min_BQ, cpus=threads, force=force,
                                                       scale_func=cov_scale_func, read_length=read_length)

    return (norm_kmer_file, norm_codon_file, contig_length_file, depth_file)


def combine_feature_tables(feature_file_list, output_dir, prefix, force=False):

    # output file handling
    if not folder_exists(output_dir):
        create_directory(output_dir)

    output_merged_file = os.path.join(output_dir, prefix + "_merged.tsv")

    try:
        emit_file_exist_warning(filename=output_merged_file, force=force)
    except Exception as e:
        _logger.info(e)
        all_feature_dfs = []
        for feature_file in feature_file_list:
            curr_df = pd.read_csv(feature_file, sep='\t', header=0, index_col=0)
            curr_df.sort_index(inplace=True)
            all_feature_dfs.append(curr_df)

        combined_df = pd.concat(all_feature_dfs, axis=1, sort=True)
        combined_df.index.name = "Contig_ID"
        combined_df.to_csv(output_merged_file, sep="\t", header=True, index=True)

    return output_merged_file
