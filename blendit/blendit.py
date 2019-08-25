# -*- coding: utf-8 -*-
"""Main module."""

import os
import sys
import click
import logging
from .utils import add_options
from .utils import setup_logging
from .commands.profile import get_codon_frequencies_for_contigs
from .commands.profile import get_kmer_frequencies_for_contigs
from .commands.profile import parallel_calculate_contig_depth_from_bam_files
from .commands.profile import blendit_profiler
from .commands.profile import blendit_get_profiles
from .commands.profile import blendit_merge_profiles
from .commands.bin import iterative_dbscan_clustering
from .commands.bin import iterative_hdbscan_clustering
from .commands.bin import iterative_optics_clustering
from .utils.embedding import compute_embeddings
from .utils.external import run_das_tool
from .utils.common import set_loglevel

import warnings
warnings.filterwarnings(action='ignore', category=DeprecationWarning)


# logger
_logger = logging.getLogger("BlendIt")
setup_logging()


def emit_subcommand_info(subcommand, loglevel):
    set_loglevel(loglevel)
    _logger.info('invoking {0} subcommand ...'.format(subcommand))


# global shared options
shared_options = [
    click.option('-p', '--prefix', help="output prefix", default="assembly", type=str, show_default=True),
    click.option('-o', '--output_dir', help="output directory", default="./blendit_results", show_default=True),
    click.option('-f', '--force', is_flag=True, default=False, help="force to overwrite the output file"),
    click.option('-l', '--loglevel', default='debug', show_default=True,
        type=click.Choice(['critical', 'error', 'warning', 'info', 'debug'])),
    click.option('-t', '--threads', type=int, default=20, show_default=True,
        help="maximum number of threads/cpus to use when available"),
    click.version_option(version="0.1.0", prog_name="blendit", message="%(prog)s, version %(version)s")]


bin_shared_arguments = [
    click.argument('kmerfreq_file', type=str),
    click.argument('codonfreq_file', type=str),
    click.argument('depth_file', type=str),
    click.argument('contig_length_file', type=str),
    click.argument('assembly', type=str)]


bin_shared_options = [
    click.option('-x', '--min_length_x', type=int, default=2000, show_default=True,
        help="minimum contig length threshold x"),
    click.option('-y', '--min_length_y', type=int, default=10000, show_default=True,
        help="minimum contig length threshold y"),
    click.option('-s', '--length_step', type=int, default=1000, show_default=True,
        help="minimum contig length increasement step"),
    click.option('-d', '--dimred', default='both', show_default=True,
        type=click.Choice(['tsne', 'umap', 'both']),
        help="dimension reduction methods, can be 'tsne', 'umap' or 'both'"),
    click.option('--dimensions', type=int, default=3, show_default=True,
        help="number of dimensions to keep for embedding"),
    click.option('--components', type=int, default=100, show_default=True,
        help="maximum PCA components to keep")]

pipe_shared_arguments = [
    click.argument('assembly', type=str),
    click.argument('bam_files', type=str, nargs=-1)]

pipe_shared_options = [
    click.option('-k', '--kmer_size', type=click.Choice(['4', '5']), default='5',
        show_default=True, help="k-mer size"),
    click.option('--kmerfreq_scale_func', type=click.Choice(['none', 'sqrt', 'cbrt', 'log10']),
        default='cbrt', show_default=True, help="k-mer freq scale function"),
    click.option('--use_kmercount', is_flag=True, default=False,
        help="use kmer count instead of kmer frequency (contig-wise/row-wise normalization)"),
    click.option('--codonfreq_scale_func', type=click.Choice(['none', 'sqrt', 'cbrt', 'log10']),
        default='cbrt', show_default=True, help="codon freq scale function"),
    click.option('--cov_scale_func', type=click.Choice(['none', 'sqrt', 'cbrt', 'log10']),
        default='log10', show_default=True, help="coverage scale function"),
    click.option('-l', '--read_length', type=int, default=250, show_default=True,
        help="read length for log-scaled transformation")]
pipe_shared_options.extend(bin_shared_options)

# entry point
@click.group()
def main(*args, **kwargs):
    """

    BlendIt: Binning metagenomic contigs via length-dependent Iterative clustering and integration

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
        -> blendit bin dbscan:     run dbscan binning
        -> blendit bin optics:     run optics binning

    \b
    -> blendit pipe COMMAND [ARGS] ...
        -> blendit pipe ph:        run all feature profilings and hdbscan clustering
        -> blendit pipe pd:        run all feature profilings and dbscan clustering

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


# pipe group
@main.group(invoke_without_command=False)
@click.pass_context
def pipe(ctx, *args, **kwargs):
    pass


# profile kmer
@profile.command()
@click.argument('assembly', type=str)
@click.option('-k', '--kmer_size', type=click.Choice(['4', '5']), default='5', show_default=True, help="k-mer size")
@click.option('--kmerfreq_scale_func', type=click.Choice(['none', 'sqrt', 'cbrt', 'log10']),
              default='cbrt', show_default=True, help="k-mer freq scale function")
@add_options(shared_options)
def kmer(assembly, kmer_size, kmerfreq_scale_func, prefix, output_dir, force, threads, loglevel, *args, **kwargs):
    """k-mer frequency profiling"""
    emit_subcommand_info("profile kmer", loglevel)

    output_dir = os.path.join(output_dir, "kmer")
    norm_kmer_file = get_kmer_frequencies_for_contigs(input_contig_file=assembly, output_dir=output_dir, prefix=prefix,
                                                      k=kmer_size, cpus=threads, force=force,
                                                      scale_func=kmerfreq_scale_func)

    return norm_kmer_file


# profile codon
@profile.command()
@click.argument('assembly', type=str)
@click.option('-g', '--genetic_code', type=int, default=11, show_default=True, help="genetic code")
@click.option('--codonfreq_scale_func', type=click.Choice(['none', 'sqrt', 'cbrt', 'log10']),
              default='cbrt', show_default=True, help="codon freq scale function")
@add_options(shared_options)
def codon(assembly, codonfreq_scale_func, prefix, output_dir, genetic_code, force, threads, loglevel, *args, **kwargs):
    """codon usage profiling"""
    emit_subcommand_info("profile codon", loglevel)

    output_dir = os.path.join(output_dir, "codon")
    norm_codon_file = get_codon_frequencies_for_contigs(input_contig_file=assembly, output_dir=output_dir,
                                                        prefix=prefix, genetic_code=genetic_code, force=force,
                                                        scale_func=codonfreq_scale_func, cpus=threads)

    return norm_codon_file


# profile cov
@profile.command()
@click.argument('bam_files', type=str, nargs=-1)
@click.option('--cov_scale_func', type=click.Choice(['none', 'sqrt', 'cbrt', 'log10']),
              default='log10', show_default=True, help="coverage scale function")
@click.option('-l', '--read_length', type=int, default=250, show_default=True, help="read length for log-scaled transformation")
@add_options(shared_options)
def cov(bam_files, cov_scale_func, read_length, prefix, output_dir, threads, force, loglevel, *args, **kwargs):
    """read coverage profiling"""
    emit_subcommand_info("profile cov", loglevel)

    output_dir = os.path.join(output_dir, "cov")
    length_file, norm_depth_file = parallel_calculate_contig_depth_from_bam_files(bam_files, output_dir=output_dir,
                                                prefix=prefix, min_read_len=30, scale_func=cov_scale_func,
                                                read_length=read_length, min_MQ=0, min_BQ=0, cpus=threads, force=force)

    return length_file, norm_depth_file


# bin hdbscan
@bin.command()
@add_options(bin_shared_arguments)
@add_options(bin_shared_options)
@add_options(shared_options)
def hdbscan(kmerfreq_file, codonfreq_file, depth_file, contig_length_file, assembly, min_length_x, min_length_y, length_step,
            threads, prefix, output_dir, force, loglevel, dimred, dimensions, components, *args, **kwargs):
    """hdbscan binning"""
    emit_subcommand_info("bin hdbscan", loglevel)

    output_dir = os.path.join(output_dir, "hdbscan")

    # get feature files
    feature_files = blendit_merge_profiles(kmerfreq_file, codonfreq_file, depth_file, output_dir, prefix, force, threads,
                           dimred, dimensions, components)

    # iterative hdbcan binning
    scaffold2bin_files, bin_folders = iterative_hdbscan_clustering(embeddings=feature_files, assembly=assembly,
                                contig_length_file=contig_length_file, output_dir=output_dir, prefix=prefix,
                                force=force, min_length_x=min_length_x, min_length_y=min_length_y, length_step=length_step)

    # run DAS_Tool
    _logger.info("integrating using DAS_Tool ...")
    proteins = codonfreq_file.replace("_codonfreq_norm.tsv", ".prot")
    labels = [os.path.basename(os.path.normpath(bin_folder)) for bin_folder in bin_folders]
    run_das_tool(bins=','.join(scaffold2bin_files), contigs=assembly, labels=','.join(labels),
                 output_dir=output_dir, output_prefix=prefix+"_dastool", proteins=proteins,
                 search_engine='usearch', write_bin_evals=1, create_plots=1, write_bins=1, threads=threads,
                 score_threshold=0.4, duplicate_penalty=0.4, megabin_penalty=0.3)


# bin dbscan
@bin.command()
@add_options(bin_shared_arguments)
@add_options(bin_shared_options)
@add_options(shared_options)
def dbscan(kmerfreq_file, codonfreq_file, depth_file, contig_length_file, assembly, min_length_x, min_length_y, length_step,
            threads, prefix, output_dir, force, loglevel, dimred, dimensions, components, *args, **kwargs):
    """dbscan binning"""
    emit_subcommand_info("bin dbscan", loglevel)

    output_dir = os.path.join(output_dir, "dbscan")

    # get feature files
    feature_files = blendit_merge_profiles(kmerfreq_file, codonfreq_file, depth_file, output_dir, prefix, force, threads,
                           dimred, dimensions, components)


    # iterative dbcan binning
    scaffold2bin_files, bin_folders = iterative_dbscan_clustering(embeddings=feature_files, assembly=assembly,
                                contig_length_file=contig_length_file, output_dir=output_dir, prefix=prefix,
                                force=force, min_length_x=min_length_x, min_length_y=min_length_y, length_step=length_step)

    # run DAS_Tool
    _logger.info("integrating using DAS_Tool ...")
    proteins = codonfreq_file.replace("_codonfreq_norm.tsv", ".prot")
    labels = [os.path.basename(os.path.normpath(bin_folder)) for bin_folder in bin_folders]
    run_das_tool(bins=','.join(scaffold2bin_files), contigs=assembly, labels=','.join(labels),
                 output_dir=output_dir, output_prefix=prefix+"_dastool", proteins=proteins,
                 search_engine='usearch', write_bin_evals=1, create_plots=1, write_bins=1, threads=threads,
                 score_threshold=0.4, duplicate_penalty=0.4, megabin_penalty=0.3)

# bin optics
@bin.command()
@add_options(bin_shared_arguments)
@add_options(bin_shared_options)
@add_options(shared_options)
def optics(kmerfreq_file, codonfreq_file, depth_file, contig_length_file, assembly, min_length_x, min_length_y, length_step,
            threads, prefix, output_dir, force, loglevel, dimred, dimensions, components, *args, **kwargs):
    """optics binning"""
    emit_subcommand_info("bin optics", loglevel)

    output_dir = os.path.join(output_dir, "optics")

    # get feature files
    feature_files = blendit_merge_profiles(kmerfreq_file, codonfreq_file, depth_file, output_dir, prefix, force, threads,
                           dimred, dimensions, components)

    # iterative dbcan binning
    scaffold2bin_files, bin_folders = iterative_optics_clustering(embeddings=feature_files, assembly=assembly,
                                contig_length_file=contig_length_file, output_dir=output_dir, prefix=prefix,
                                force=force, min_length_x=min_length_x, min_length_y=min_length_y, length_step=length_step)
    #                           min_cluster_size=10)

    # run DAS_Tool
    _logger.info("integrating using DAS_Tool ...")
    proteins = codonfreq_file.replace("_codonfreq_norm.tsv", ".prot")
    labels = [os.path.basename(os.path.normpath(bin_folder)) for bin_folder in bin_folders]
    run_das_tool(bins=','.join(scaffold2bin_files), contigs=assembly, labels=','.join(labels),
                 output_dir=output_dir, output_prefix=prefix+"_dastool", proteins=proteins,
                 search_engine='usearch', write_bin_evals=1, create_plots=1, write_bins=1, threads=threads,
                 score_threshold=0.4, duplicate_penalty=0.4, megabin_penalty=0.3)


# pipe ph
@pipe.command()
@add_options(pipe_shared_arguments)
@add_options(pipe_shared_options)
@add_options(shared_options)
@click.pass_context
def ph(ctx, assembly, bam_files, prefix, output_dir, kmer_size, kmerfreq_scale_func, use_kmercount, codonfreq_scale_func, cov_scale_func,
       min_length_x, min_length_y, length_step, threads, dimred, dimensions, components, read_length, force, loglevel,
       *args, **kwargs):
    """run feature profiling and hdbscan clustering pipeline """

    emit_subcommand_info("pipe ph", loglevel)

    kmerfreq_file, codonfreq_file, contig_length_file, depth_file = blendit_get_profiles(assembly=assembly, bam_files=bam_files,
                     prefix=prefix, output_dir=output_dir, kmer_size=kmer_size, kmerfreq_scale_func=kmerfreq_scale_func,
                     codonfreq_scale_func=codonfreq_scale_func, genetic_code=11, cov_scale_func=cov_scale_func, min_read_len=30,
                     min_MQ=0, min_BQ=0, threads=threads, read_length=read_length, force=force)

    ctx.invoke(hdbscan, kmerfreq_file=kmerfreq_file, codonfreq_file=codonfreq_file, depth_file=depth_file,
               contig_length_file=contig_length_file, assembly=assembly, min_length_x=min_length_x, min_length_y=min_length_y,
                length_step=length_step, threads=threads, prefix=prefix, output_dir=output_dir, force=force,
               loglevel=loglevel, dimred=dimred, dimensions=dimensions, components=components, *args, **kwargs)


# pipe pd
@pipe.command()
@add_options(pipe_shared_arguments)
@add_options(pipe_shared_options)
@add_options(shared_options)
@click.pass_context
def pd(ctx, assembly, bam_files, prefix, output_dir, kmer_size, kmerfreq_scale_func, use_kmercount, codonfreq_scale_func, cov_scale_func,
       min_length_x, min_length_y, length_step, threads, dimred, dimensions, components, read_length, force, loglevel,
       *args, **kwargs):
    """run feature profiling and dbscan clustering pipeline """

    emit_subcommand_info("pipe pd", loglevel)

    kmerfreq_file = ctx.invoke(kmer, assembly=assembly, kmer_size=kmer_size, kmerfreq_scale_func=kmerfreq_scale_func, use_kmercount=use_kmercount,
                               prefix=prefix, output_dir=output_dir, force=force, cpus=threads, loglevel=loglevel,
                               *args, **kwargs)
    #kmerfreq_file = os.path.join(output_dir, "kmer", prefix + "_kmer_norm.tsv")
    codonfreq_file = ctx.invoke(codon, assembly=assembly, codonfreq_scale_func=codonfreq_scale_func, prefix=prefix, output_dir=output_dir,
                                force=force, loglevel=loglevel, *args, **kwargs)
    #codonfreq_file = os.path.join(output_dir, "codon", prefix + "_codon_norm.tsv")
    contig_length_file, depth_file = ctx.invoke(cov, bam_files=bam_files, cov_scale_func=cov_scale_func, prefix=prefix, output_dir=output_dir,
                                                threads=threads, force=force, loglevel=loglevel, *args, **kwargs)
    #depth_file = os.path.join(output_dir, "cov", prefix + "_cov_norm.tsv")
    contig_length_file = os.path.join(output_dir, "cov", prefix + "_cov_contig_length.tsv")

    ctx.invoke(dbscan, kmerfreq_file=kmerfreq_file, codonfreq_file=codonfreq_file, depth_file=depth_file,
               contig_length_file=contig_length_file, assembly=assembly, min_length_x=min_length_x, min_length_y=min_length_y,
                length_step=length_step, threads=threads, prefix=prefix, output_dir=output_dir, force=force,
               loglevel=loglevel, dimred=dimred, dimensions=dimensions, components=components, *args, **kwargs)


if __name__ == "__main__":
    sys.exit(main())
