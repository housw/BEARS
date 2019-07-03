# -*- coding: utf-8 -*-
"""Main module."""

import os
import sys
import click
import logging
import subprocess
import numpy as np
import pandas as pd
from .utils import add_options
from .utils import setup_logging
from .utils import CommandWrapper
from .utils import SpecialHelpOrder
from .utils import run_prodigal
from .utils import get_prefix
from .utils import command_logger
from .utils.codon import get_codon_frequencies_for_contigs
from .utils.common import normalizer
from .utils.common import combine_feature_tables
from .utils.coordinate import compute_PCA_tSNE_coordinates
from .utils.coordinate import compute_tSNE_coordinates
from .utils.kmer import get_kmer_counts_for_contigs
from .utils.cov import calculate_contig_depth_from_bam_files
from .utils.cov import normalize_contig_depth
from .utils.hdbscan import hdbscan_clustering
from .utils.hdbscan import generate_bins


_logger = logging.getLogger("bears")


def emit_subcommand_info(subcommand, loglevel):
    setup_logging(loglevel)
    _logger.info('invoking {0} subcommand'.format(subcommand))


# entry point
@click.group()
def main(*args, **kwargs):
    """BEARS: Binning mEtagenomic contigs via iterAtive hieRarchical dbScan. 

    The program implements the `profile` and `bin` steps of 
    metagenomic data analysis. You can invoke each step using the 
    following subcommands: 

    \b
    -> bears profile COMMAND [ARGS] ...
    \b
      -> bears profile codon:  run codon usage profiling
      -> bears profile cov:    run read coverage profiling
      -> bears profile kmer:   run k-mer frequency profiling
 
    \b
    -> bears bin COMMAND [ARGS] ...
      -> bears bin hdbscan:    run hdbscan binning
      -> bears bin dbscanpp:   run dbscan++ binning

    """
    pass


# profile group
@main.group(invoke_without_command=False)
@click.pass_context
def profile(ctx, *args, **kwargs):
     pass


# bin group
@main.group(invoke_without_command=False)
@click.pass_context
def bin(ctx, *args, **kwargs):
    pass


# global shared options
shared_options = [
    click.option('-p', '--prefix', help="output prefix", type=str), 
    click.option('-o', '--output_dir', help="output directory", default="./bears_results", show_default=True),
    click.option('-f', '--force', is_flag=True, default=False, help="force to overwrite the output file"), 
    click.option('-l', '--loglevel', default='debug', show_default=True,
        type=click.Choice(['critical', 'error', 'warning', 'info', 'debug'])),
    click.version_option(version="0.1.0", prog_name="bears", message="%(prog)s, version %(version)s")
]


# profile kmer
@profile.command()
@click.argument('assembly', type=str)
@click.option('-k', '--kmer_size', type=click.Choice(['4', '5']), default='5', 
    show_default=True, help="k-mer size")
@click.option('-c', '--cpus', type=int, default=20, show_default=True, help="number of cores to use for kmer counting")
@add_options(shared_options)
def kmer(assembly, kmer_size, prefix, output_dir, force, cpus, loglevel, *args, **kwargs):
    """k-mer frequency profiling"""
    emit_subcommand_info("profile kmer", loglevel)

    if not prefix:
        prefix = get_prefix(assembly)

    # run kmer freq calculation
    kmer_freq_file = get_kmer_counts_for_contigs(input_contig_file=assembly, output_dir=output_dir, prefix=prefix,
                                                 k=kmer_size, cpus=cpus, force=force)

    # normalization
    norm_kmer_file = normalizer(input_freq_file=kmer_freq_file, output_dir=output_dir, prefix=prefix+"_kmerfreq", scale_func=np.cbrt)

    # t-SNE
    #tSNE_freq_file = compute_tSNE_coordinates(freq_file=norm_kmer_file, output_dir=output_dir,
    #                                          prefix=prefix+"_kmerfreq_norm", threads=cpus, force=force,
    #                                          columns=['KmerFreq_tSNE_X', 'KmerFreq_tSNE_Y'])



# profile codon
@profile.command()
@click.argument('assembly', type=str)
@add_options(shared_options)
def codon(assembly, prefix, output_dir, force, loglevel, *args, **kwargs):
    """codon usage profiling"""
    emit_subcommand_info("profile codon", loglevel)

    if not prefix:
        prefix = get_prefix(assembly)

    # run prodigal
    prodigal_file, prodigal_gene, prodigal_prot = \
        run_prodigal(assembly, prefix=prefix, output_dir=output_dir, output_fmt='gbk', flags=['m'], force=force)

    # run codon usage calculation
    codon_freq_file = get_codon_frequencies_for_contigs(input_prodigal_nucl_file=prodigal_gene, output_dir=output_dir,
                                                        prefix = prefix, genetic_code=11, force=force)

    # normalization
    norm_codon_file = normalizer(input_freq_file=norm_codon_file, output_dir=output_dir, prefix=prefix+"_codonfreq", scale_func=np.cbrt)

    # t-SNE
    #tSNE_freq_file = compute_tSNE_coordinates(freq_file=norm_freq_file, output_dir=output_dir,
    #                                          prefix=prefix+"_codonfreq_norm", threads=20, force=force,
    #                                          columns=['CodonFreq_tSNE_X', 'CodonFreq_tSNE_Y'])

# profile cov
@profile.command()
@click.argument('bam_files', type=str, nargs=-1)
@click.option('-l', '--read_length', type=int, default=250, show_default=True, help="read length for log-scaled transformation")
@add_options(shared_options)
def cov(bam_files, prefix, output_dir, force, loglevel, *args, **kwargs):
    """read coverage profiling"""
    emit_subcommand_info("profile cov", loglevel)

    if not prefix:
        prefix = "bamcov"

    # run bamcov
    length_file, depth_file = calculate_contig_depth_from_bam_files(bam_files, output_dir=output_dir,
                                                                    output_prefix=prefix, min_read_len=30,
                                                                    min_MQ=0, min_BQ=0, force=force)
    # normalization
    norm_depth_file = normalize_contig_depth(length_file, depth_file, output_dir, prefix, scale_func=np.log10, read_length=250)

    # t-SNE
    #tSNE_freq_file = compute_tSNE_coordinates(freq_file=norm_depth_file, output_dir=output_dir,
    #                                          prefix=prefix+"_depth_norm", threads=20, force=force,
    #                                          columns=['Depth_tSNE_X', 'Depth_tSNE_Y'])

# bin hdbscan
@bin.command()
@click.argument('kmerfreq_file', type=str)
@click.argument('codonfreq_file', type=str)
@click.argument('depth_file', type=str)
@click.argument('contig_length_file', type=str)
@click.argument('assembly', type=str)
@click.option('-x', '--min_length_x', type=int, default=2000, show_default=True,
    help="minimum contig length threshold x")
@click.option('-y', '--min_length_y', type=int, default=10000, show_default=True, 
    help="minimum contig length threshold y")
@click.option('-s', '--length_step', type=int, default=1000, show_default=True,
    help="minimum contig length increasement step")
@add_options(shared_options)
def hdbscan(kmerfreq_file, codonfreq_file, depth_file, contig_length_file, assembly, min_length_x, min_length_y, length_step, prefix, output_dir, force, loglevel, *args, **kwargs):
    """hdbscan binning"""
    emit_subcommand_info("bin hdbscan", loglevel)

    if not prefix:
        prefix = "hdbscan"

    # merge kmer, codon and depth profiles
    merged_profile = combine_feature_tables([kmerfreq_file, codonfreq_file, depth_file], output_dir=output_dir, prefix=prefix, force=force)

    # run t-SNE
    tsne_file = compute_PCA_tSNE_coordinates(merged_profile, output_dir, prefix+"_merged", threads=20, n_components=100)

    # run hdbscan
    length_df = pd.read_csv(contig_length_file, sep="\t", index_col=0, header=0)
    for min_length in range(min_length_x, min_length_y+1, length_step):
        print(min_length)
        drop_length_df = length_df[length_df.Length < min_length]
        drop_contig_list = drop_length_df.index
        hdbscan_cluster_file = hdbscan_clustering(tsne_file, output_dir, prefix+"_merged_tSNE_min_{}".format(min_length), excluding_contig_file_or_list=drop_contig_list, min_cluster_size=10)
        # generate bins
        generate_bins(hdbscan_cluster_file, assembly, output_dir, bin_folder=prefix+"_merged_tSNE_min_{}".format(min_length))



# https://arxiv.org/pdf/1810.13105.pdf
# bin dbscanpp
@bin.command()
@click.argument('tsne_profile', type=str)
@click.option('-x', '--min_length_x', type=int, default=2000, show_default=True,
    help="minimum contig length threshold x")
@click.option('-y', '--min_length_y', type=int, default=10000, show_default=True, 
    help="minimum contig length threshold y")
@add_options(shared_options)
def dbscanpp(tsne_profile, min_length_x, min_length_y, prefix, output_dir, length_threshold, force, loglevel, *args, **kwargs):
    """dbscan++ binning"""
    emit_subcommand_info("bin dbscanpp", loglevel)
    pass


if __name__ == "__main__":
    sys.exit(main()) 
