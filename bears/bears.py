# -*- coding: utf-8 -*-
"""Main module."""

import os
import sys
import click
import logging
import subprocess
import numpy as np
from .utils import add_options
from .utils import setup_logging
from .utils import CommandWrapper
from .utils import SpecialHelpOrder
from .utils import run_prodigal
from .utils import get_prefix
from .utils import command_logger
from .utils.codon import get_codon_frequencies_for_contigs
from .utils.common import normalizer
from .utils.coordinate import compute_tSNE_coordinates
from .utils.kmer import get_kmer_counts_for_contigs

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
      -> bears bin umap:       run umap binning

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
    click.option('-l', '--length_threshold', type=int, default=2000, show_default=True, 
    help="minimum contig length threshold"),
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
def kmer(assembly, kmer_size, prefix, output_dir, length_threshold, force, cpus, loglevel, *args, **kwargs):
    """k-mer frequency profiling"""
    emit_subcommand_info("profile kmer", loglevel)

    if not prefix:
        prefix = get_prefix(assembly)

    # run kmer freq calculation
    kmer_freq_file = get_kmer_counts_for_contigs(input_contig_file=assembly, output_dir=output_dir, prefix=prefix,
                                                 k=kmer_size, cpus=cpus, force=False)

    # normalization
    norm_freq_file = normalizer(input_freq_file=kmer_freq_file, output_dir=output_dir, prefix=prefix+"_kmerfreq", scale_func=np.cbrt)

    # t-SNE
    tSNE_freq_file = compute_tSNE_coordinates(codon_freq_file=norm_freq_file, output_dir=output_dir,
                                              prefix=prefix+"_kmerfreq_norm", threads=20, force=force,
                                              columns=['KmerFreq_tSNE_X', 'KmerFreq_tSNE_Y'])



# profile codon
@profile.command()
@click.argument('assembly', type=str)
@add_options(shared_options)
def codon(assembly, prefix, output_dir, length_threshold, force, loglevel, *args, **kwargs):
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
    norm_freq_file = normalizer(input_freq_file=codon_freq_file, output_dir=output_dir, prefix=prefix+"_codonfreq", scale_func=np.cbrt)

    # t-SNE
    tSNE_freq_file = compute_tSNE_coordinates(codon_freq_file=norm_freq_file, output_dir=output_dir,
                                              prefix=prefix+"_codonfreq_norm", threads=20, force=force,
                                              columns=['CodonFreq_tSNE_X', 'CodonFreq_tSNE_Y'])

# profile cov
@profile.command()
@click.argument('bam_files', type=str, nargs=-1)
@add_options(shared_options)
def cov(bam_files, prefix, output_dir, length_threshold, force, loglevel, *args, **kwargs):
    """read coverage profiling"""
    emit_subcommand_info("profile cov", loglevel)
    pass


# bin hdbscan
@bin.command()
@click.argument('tsne_profile', type=str)
@click.option('-x', '--min_length_x', type=int, default=2000, show_default=True,
    help="minimum contig length threshold x")
@click.option('-y', '--min_length_y', type=int, default=10000, show_default=True, 
    help="minimum contig length threshold y")
@add_options(shared_options)
def hdbscan(tsne_profile, min_length_x, min_length_y, prefix, output_dir, length_threshold, force, loglevel, *args, **kwargs):
    """hdbscan binning"""
    emit_subcommand_info("bin hdbscan", loglevel)
    pass


# bin umap
@bin.command()
@click.argument('tsne_profile', type=str)
@click.option('-x', '--min_length_x', type=int, default=2000, show_default=True,
    help="minimum contig length threshold x")
@click.option('-y', '--min_length_y', type=int, default=10000, show_default=True, 
    help="minimum contig length threshold y")
@add_options(shared_options)
def umap(tsne_profile, min_length_x, min_length_y, prefix, output_dir, length_threshold, force, loglevel, *args, **kwargs):
    """umap binning"""
    emit_subcommand_info("bin umap", loglevel)
    pass


if __name__ == "__main__":
    sys.exit(main()) 
