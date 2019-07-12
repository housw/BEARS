# -*- coding: utf-8 -*-

import os
import logging
import subprocess
from .common import command_logger
from .common import CommandException
from .common import CommandWrapper
from .common import folder_exists
from .common import create_directory


_logger = logging.getLogger("bears")


def run_prodigal(assembly, prefix, output_dir, output_fmt='gbk', flags = ['m'], force=False):
    """
    :param assembly:
    :param prefix:
    :param output_dir:
    :param output_fmt:
    :return:
    """

    # output file handling
    if not folder_exists(output_dir):
        create_directory(output_dir)
    output_file = os.path.join(output_dir, prefix + "." + output_fmt)
    output_gene = os.path.join(output_dir, prefix + ".gene")
    output_prot = os.path.join(output_dir, prefix + ".prot")
    if os.path.exists(output_file):
        if not force:
            err_msg = "{0} exists, use --force if you want to re-run prodigal".format(output_file)
            _logger.error(err_msg)
            raise CommandException(err_msg)
        else:
            warn_msg = "re-run prodigal, {0} will be over-writen".format(output_file)
            _logger.warn(warn_msg)


    # run prodigal
    prodigal = CommandWrapper(name="prodigal",
                              arguments=[],
                              options={'i': assembly, 'o': output_file, 'd': output_gene, 'a': output_prot,
                                       'f': output_fmt},
                              flags=flags)
    prodigal.construct_command(option_prefix_char="-", flag_prefix_char="-")
    prodigal.run()

    return output_file, output_gene, output_prot


def run_bamcov(bam_file, depth_prefix, output_dir, min_read_len=30, min_MQ=0, min_BQ=0, force=False, flags=[]):

    # output file handling
    if not folder_exists(output_dir):
        create_directory(output_dir)
    output_depth = os.path.join(output_dir, depth_prefix + ".bamcov")

    if os.path.exists(output_depth):
        if not force:
            err_msg = "{0} exists, use --force if you want to re-run bamcov".format(output_depth)
            _logger.error(err_msg)
            raise CommandException(err_msg)
        else:
            warn_msg = "re-run bamcov, {0} will be over-writen".format(output_depth)
            _logger.warn(warn_msg)

    # run bamcov
    bamcov = CommandWrapper(name="bamcov",
                              arguments=[bam_file],
                              options={'output': output_depth, 'min-read-len': str(min_read_len),
                                       'min-MQ': str(min_MQ), 'min-BQ': str(min_BQ)},
                              flags=flags)
    bamcov.construct_command(option_prefix_char="--", flag_prefix_char="-")
    bamcov.run()

    return output_depth


def run_das_tool(bins, contigs, labels, output_dir, output_prefix, proteins=None, search_engine='usearch',
                 write_bin_evals=1, create_plots=1, write_bins=1, threads=1, db_directory=None,
                 score_threshold=0.5, duplicate_penalty=0.6, megabin_penalty=0.5):

    # output file handling
    if not folder_exists(output_dir):
        create_directory(output_dir)
    outputbasename = os.path.join(output_dir, output_prefix)

    options={'bins':bins, 'contigs':contigs, 'outputbasename':outputbasename,
             'labels':labels, 'search_engine':search_engine, 'write_bin_evals':write_bin_evals, 
             'create_plots':create_plots, 'write_bins':write_bins, 'threads':threads,
             'score_threshold':score_threshold, 'duplicate_penalty':duplicate_penalty, 'megabin_penalty':megabin_penalty}
    if proteins:
        options['proteins'] = proteins
    if db_directory:
        options['db_directory'] = db_directory

    # run DAS_Tool
    dastool = CommandWrapper(name="DAS_Tool", options=options) 
    dastool.construct_command(option_prefix_char="--")
    dastool.run()

