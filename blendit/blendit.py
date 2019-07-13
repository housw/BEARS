# -*- coding: utf-8 -*-
"""Main module."""

import os
import sys
import click
import logging
import numpy as np
import pandas as pd
from .utils import add_options
from .utils import setup_logging
from .utils import run_prodigal
from .utils import get_prefix
from .commands.profile import get_codon_frequencies_for_contigs
from .commands.profile import get_kmer_counts_for_contigs
from .commands.profile import calculate_contig_depth_from_bam_files
from .commands.profile import normalize_contig_depth
from .commands.bin import hdbscan_clustering
from .commands.bin import generate_bins
from .utils.common import normalizer
from .utils.common import combine_feature_tables
from .utils.coordinate import compute_PCA_tSNE_coordinates
from .utils.coordinate import compute_PCA_UMAP_coordinates
from .utils.common import scaffolds_to_bins
from .utils.external import run_das_tool


_logger = logging.getLogger("blendit")


def emit_subcommand_info(subcommand, loglevel):
    setup_logging(loglevel)
    _logger.info('invoking {0} subcommand ...'.format(subcommand))


# entry point
@click.group()
def main(*args, **kwargs):
    """Blendit: Binning metagenomic contigs via length-dependent iterative clustering and integration

    The program implements the `profile` and `bin` steps of 
    metagenomic data analysis. You can invoke each step using the 
    following subcommands: 

    \b
    -> blendit profile COMMAND [ARGS] ...
    \b
      -> blendit profile codon:  run codon usage profiling
      -> blendit profile cov:    run read coverage profiling
      -> blendit profile kmer:   run k-mer frequency profiling
 
    \b
    -> blendit bin COMMAND [ARGS] ...
      -> blendit bin hdbscan:    run hdbscan binning

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
    click.option('-p', '--prefix', help="output prefix", default="assembly", type=str, show_default=True),
    click.option('-o', '--output_dir', help="output directory", default="./blendit_results", show_default=True),
    click.option('-f', '--force', is_flag=True, default=False, help="force to overwrite the output file"), 
    click.option('-l', '--loglevel', default='debug', show_default=True,
        type=click.Choice(['critical', 'error', 'warning', 'info', 'debug'])),
    click.version_option(version="0.1.0", prog_name="blendit", message="%(prog)s, version %(version)s")
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

    prefix = prefix + "_kmer"
    output_dir = os.path.join(output_dir, "kmer")

    # run kmer freq calculation
    _logger.info("calculating kmer frequency ...")
    kmer_freq_file = get_kmer_counts_for_contigs(input_contig_file=assembly, output_dir=output_dir, prefix=prefix,
                                                 k=kmer_size, cpus=cpus, force=force)

    # normalization
    _logger.info("normalizing kmer frequency ...")
    norm_kmer_file = normalizer(input_freq_file=kmer_freq_file, output_dir=output_dir, prefix=prefix, scale_func=np.cbrt)


# profile codon
@profile.command()
@click.argument('assembly', type=str)
@add_options(shared_options)
def codon(assembly, prefix, output_dir, force, loglevel, *args, **kwargs):
    """codon usage profiling"""
    emit_subcommand_info("profile codon", loglevel)

    prefix = prefix + "_codon"
    output_dir = os.path.join(output_dir, "codon")


    # run prodigal
    _logger.info("predicting coding genes using prodigal ...")
    prodigal_file, prodigal_gene, prodigal_prot = \
        run_prodigal(assembly, prefix=prefix, output_dir=output_dir, output_fmt='gbk', flags=['m'], force=force)

    # run codon usage calculation
    _logger.info("calculating codon frequency ...")
    codon_freq_file = get_codon_frequencies_for_contigs(input_prodigal_nucl_file=prodigal_gene, output_dir=output_dir,
                                                        prefix = prefix, genetic_code=11, force=force)

    # normalization
    _logger.info("normalizing codon frequency ...")
    norm_codon_file = normalizer(input_freq_file=codon_freq_file, output_dir=output_dir, prefix=prefix, scale_func=np.cbrt)


# profile cov
@profile.command()
@click.argument('bam_files', type=str, nargs=-1)
@click.option('-l', '--read_length', type=int, default=250, show_default=True, help="read length for log-scaled transformation")
@add_options(shared_options)
def cov(bam_files, prefix, output_dir, force, loglevel, *args, **kwargs):
    """read coverage profiling"""
    emit_subcommand_info("profile cov", loglevel)

    prefix = prefix + "_cov"
    output_dir = os.path.join(output_dir, "cov")

    # run bamcov
    _logger.info("calculating contig coverage using bamcov ...")
    length_file, depth_file = calculate_contig_depth_from_bam_files(bam_files, output_dir=output_dir,
                                                                    output_prefix=prefix, min_read_len=30,
                                                                    min_MQ=0, min_BQ=0, force=force)

    # normalization
    _logger.info("normalizing contig coverage ...")
    norm_depth_file = normalize_contig_depth(length_file, depth_file, output_dir, prefix, scale_func=np.log10, read_length=250)


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
@click.option('-t', '--threads', type=int, default=20, show_default=True,
    help="maximum number of threads to use when available")
@click.option('-d', '--dimred', default='both', show_default=True,
             type=click.Choice(['tsne', 'umap', 'both']),
             help="dimension reduction methods, can be 'tsne', 'umap' or 'both'")
@click.option('--dimensions', type=int, default=3, show_default=True,
             help="number of dimensions to keep for embedding")
@click.option('--components', type=int, default=100, show_default=True,
             help="maximum PCA components to keep")
@add_options(shared_options)
def hdbscan(kmerfreq_file, codonfreq_file, depth_file, contig_length_file, assembly, min_length_x, min_length_y, length_step,
            threads, prefix, output_dir, force, loglevel, dimred, dimensions, components, *args, **kwargs):
    """hdbscan binning"""
    emit_subcommand_info("bin hdbscan", loglevel)

    prefix = prefix + "_hdbscan"
    output_dir = os.path.join(output_dir, "hdbscan")


    # merge kmer, codon and depth profiles
    _logger.info("combining kmer, codon and coverage profiles ...")
    merged_profile = combine_feature_tables([kmerfreq_file, codonfreq_file, depth_file], output_dir=output_dir, prefix=prefix, force=force)
    #merged_profile = os.path.join(output_dir, prefix+"_merged.tsv")

    feature_files = []
    if dimred in ('tsne', 'both'):
        # run t-SNE
        _logger.info("computing t-SNE coordinates after PCA ...")
        tsne_file = compute_PCA_tSNE_coordinates(merged_profile, output_dir, prefix+"_merged_{}d".format(dimensions), threads=threads, n_components=dimensions, pca_components=components)
        #tsne_file =  os.path.join(output_dir, prefix+"_merged_{}d.tsne".format(dimensions))
        feature_files.append(tsne_file)
    if dimred in ('umap', 'both'):
        # run UMAP
        _logger.info("computing UMAP coordinates after PCA ...")
        umap_file = compute_PCA_UMAP_coordinates(merged_profile, output_dir, prefix+"_merged_{}d".format(dimensions), n_components=dimensions, pca_components=components)
        #umap_file = os.path.join(output_dir, prefix + "_merged_{}d.umap".format(dimensions))
        feature_files.append(umap_file)


    scaffold2bin_files = []
    bin_folders = []
    for feature_file in feature_files:
        dimred_method = feature_file.split('.')[-1]
        _logger.info("clustering using hdbscan based on {0} feature file ...".format(dimred_method))
        length_df = pd.read_csv(contig_length_file, sep="\t", index_col=0, header=0)
        for min_length in range(min_length_x, min_length_y+1, length_step):
            _logger.info("clustering using hdbscan with contigs >= {l} ...".format(l=min_length))
            drop_length_df = length_df[length_df.Length < min_length]
            drop_contig_list = drop_length_df.index
            hdbscan_cluster_file = hdbscan_clustering(feature_file, output_dir, prefix+"_merged_{0}_min_{1}".format(dimred_method, min_length), excluding_contig_file_or_list=drop_contig_list, min_cluster_size=10)
            # generate bins
            _logger.info("generating bins for contigs >= {l} ...".format(l=min_length))
            generate_bins(hdbscan_cluster_file, assembly, output_dir, bin_folder=prefix+"_merged_{0}_min_{1}".format(dimred_method, min_length))
            # generate scaffold2bin files
            _logger.info("generating scaffold2bin file for contigs >= {0} ...".format(min_length))
            bin_folder = os.path.join(output_dir, prefix+"_merged_{0}_min_{1}".format(dimred_method, min_length))
            scaffold2bin = bin_folder + "_scaffold2bin.tsv"
            bin_folders.append(bin_folder)
            scaffold2bin_files.append(scaffold2bin)
            scaffolds_to_bins(bin_folder, scaffold2bin, suffix="fa")
        

    # run DAS_Tool
    _logger.info("integrating using DAS_Tool ...")
    proteins = codonfreq_file.replace("_norm.tsv", ".prot")
    labels = [os.path.basename(os.path.normpath(bin_folder)) for bin_folder in bin_folders]
    run_das_tool(bins=','.join(scaffold2bin_files), contigs=assembly, labels=','.join(labels),
                 output_dir=output_dir, output_prefix=prefix+"_dastool", proteins=proteins,
                 search_engine='usearch', write_bin_evals=1, create_plots=1, write_bins=1, threads=threads,
                 score_threshold=0.4, duplicate_penalty=0.4, megabin_penalty=0.3)

if __name__ == "__main__":
    sys.exit(main()) 
