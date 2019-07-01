# -*- coding: utf-8 -*-

"""Console script for bears."""

import os
import sys
import click
import logging
from .bears import * 
from .utils import add_options
from .utils import setup_logging
from .utils import SpecialHelpOrder

_logger = logging.getLogger(__name__)
def emit_subcommand_info(subcommand, loglevel):
    setup_logging(loglevel)
    _logger.info('invoking {0} subcommand'.format(subcommand))


# entry point
@click.group()
def main(**kwargs):
    """BEARS: Binning mEtagenomic contigs via iterAtive hieRarchical dbScan. 

    The program implements the `profile` and `bin` steps of 
    metagenomic data analysis. You can invoke each step using the 
    following subcommands: 

    \b
    > bears profile COMMAND [ARGS] ...
    \b
      - bears profile codon:  run codon usage profiling
      - bears profile cov:    run read coverage profiling
      - bears profile kmer:   run k-mer frequency profiling
 
    \b
    > bears bin COMMAND [ARGS] ...
      - bears bin hdbscan:    run hdbscan binning
      - bears bin umap:       run umap binning

    """
    pass


# profile group
@main.group(invoke_without_command=False)
@click.pass_context
def profile(ctx, **kwargs):
     pass


# bin group
@main.group(invoke_without_command=False)
@click.pass_context
def bin(ctx, **kwargs):
    pass


# global shared options
shared_options = [
    click.option('-p', '--prefix', help="output prefix", type=str), 
    click.option('-o', '--output_dir', help="output directory", default="./", show_default=True), 
    click.option('-f', '--force', is_flag=True, default=False, help="force to overwrite the output file"), 
    click.option('-l', '--length_threshold', type=int, default=2000, show_default=True, 
    help="minimum contig length threshold"),
    click.option('-l', '--loglevel', default='info', show_default=True,
        type=click.Choice(['critical', 'error', 'warning', 'info', 'debug'])),
    click.version_option(version="0.1.0", prog_name="bears", message="%(prog)s, version %(version)s")
]


# profile kmer
@profile.command()
@click.argument('assembly', type=str)
@click.option('-k', '--kmer_size', type=click.Choice(['4', '5']), default='5', 
    show_default=True, help="k-mer size")
@add_options(shared_options)
def kmer(assembly, kmer_size, prefix, output_dir, length_threshold, force, loglevel):
    """k-mer frequency profiling"""
    emit_subcommand_info("profile kmer", loglevel)
    pass

# profile codon
@profile.command()
@click.argument('assembly', type=str)
@add_options(shared_options)
def codon(assembly, prefix, output_dir, length_threshold, force, loglevel):
    """codon usage profiling"""
    emit_subcommand_info("profile codon", loglevel)
    pass


# profile cov
@profile.command()
@click.argument('bam_files', type=str, nargs=-1)
@add_options(shared_options)
def cov(bam_files, prefix, output_dir, length_threshold, force, loglevel):
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
def hdbscan(tsne_profile, min_length_x, min_length_y, prefix, output_dir, length_threshold, force, loglevel):
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
def umap(tsne_profile, min_length_x, min_length_y, prefix, output_dir, length_threshold, force, loglevel):
    """umap binning"""
    emit_subcommand_info("bin umap", loglevel)
    pass


if __name__ == "__main__":
    sys.exit(main()) 
