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


shared_options = [
    click.option('-p', '--prefix', help="output prefix", type=str), 
    click.option('-o', '--output_dir', help="output directory", default="./", show_default=True), 
    click.option('-f', '--force', is_flag=True, default=False, help="force to overwrite the output file"), 
    click.option('-l', '--loglevel', default='info', show_default=True,
        type=click.Choice(['critical', 'error', 'warning', 'info', 'debug'])),
    click.version_option(version="0.1.0", prog_name="bears", message="%(prog)s, version %(version)s")
]


@click.group(cls=SpecialHelpOrder)
def main(**kwargs):
    pass


@main.command(help_priority=1)
@click.argument('assembly', type=str)
@click.argument('bam_files', type=str, nargs=-1)
@click.option('-k', '--kmer_size', type=click.Choice(['4', '5']), default='5', 
    show_default=True, help="k-mer size")
@click.option('-l', '--length_threshold', type=int, default=2000, show_default=True, 
    help="minimum contig length threshold")
@add_options(shared_options)
def profile(assembly, bam_files, prefix, output_dir, force, loglevel):
    emit_subcommand_info("profile", loglevel)
    #output_file = make_output_file(input_file, prefix, output_dir, force, suffix=".txt")
    # kmer profile 
    pass
    # coverage profile
    pass


@main.command(help_priority=2)
@click.argument('tSNE_profile', type=str)
@click.option('-x', '--min_length_x', type=int, default=2000, show_default=True,
    help="minimum contig length threshold x")
@click.option('-y', '--min_length_y', type=int, default=10000, show_default=True, 
    help="minimum contig length threshold y")
@add_options(shared_options)
def bin(tSNE_profile, min_length_x, min_length_y, prefix, output_dir, force, loglevel):
    emit_subcommand_info("profile", loglevel)


if __name__ == "__main__":
    sys.exit(main()) 
