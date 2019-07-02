# -*- coding: utf-8 -*-


import os
import sys
import click
import errno
import logging
import subprocess
import numpy as np
import pandas as pd
from functools import wraps
from sklearn import preprocessing


_logger = logging.getLogger("bears")


# credit: https://github.com/pallets/click/issues/108
def add_options(options):
    def _add_options(func):
        for option in reversed(options):
            func = option(func)
        return func
    return _add_options


def make_output_file(input_file, prefix=None, output_dir="./", force=False, suffix=".txt"):
    """make output_file, check existence"""

    # input and output handeling
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not prefix:
        basename = os.path.basename(os.path.normpath(input_file))
        _logger.debug("output basename is {}".format(basename))
        if "." in basename:
            prefix, ext = os.path.splitext(basename)
        else:
            prefix, ext = basename, "None"
        # e.g., remove .fq.gz
        if ext in (".gz", ".gzip", ".bz2", ".bzip2", ".zip"):
            prefix = os.path.splitext(prefix)[0]
    _logger.info("output prefix is {}".format(prefix))
    out_file = os.path.join(output_dir, prefix + suffix)
    _logger.info("output file is {}".format(out_file))
    if os.path.exists(out_file):
        if force:
            _logger.warning("output file exists, will be overwritten!")
        else:
            err_msg = "output file detected, please backup it at first!\n\n"
            _logger.error(err_msg)
            raise click.UsageError(message=err_msg)
    return out_file


def setup_logging(loglevel):
    """Setup basic loggings
    Args:
      loglevel (str): minimum loglevel for emitting messages
    """

    loglevel = {
        'critical': logging.CRITICAL,
        'error'   : logging.ERROR,
        'warning' : logging.WARNING,
        'info'    : logging.INFO,
        'debug'   : logging.DEBUG,
    }.get(loglevel, logging.DEBUG)

    logformat = "[%(asctime)s] [%(levelname)s] %(name)s:%(message)s"

    logging.basicConfig(level=loglevel, format=logformat, datefmt="%Y-%m-%d %H:%M:%S")


def command_logger(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        arguments = " ".join(args)
        for k, v in kwargs.items():
            if len(k) == 1:
                arguments += " -" + k + " " + str(v)
            else:
                arguments += " --" + k + " " + str(v)
        _logger.info("running external command: " + func.__name__ + " "+ arguments)
        return func(*args, **kwargs)
    return wrapper


# credit: 
# https://stackoverflow.com/questions/47972638/how-can-i-define-the-order-of-click-sub-commands-in-help/47984810#47984810
class SpecialHelpOrder(click.Group):

    def __init__(self, *args, **kwargs):
        self.help_priorities = {}
        super(SpecialHelpOrder, self).__init__(*args, **kwargs)

    def get_help(self, ctx):
        self.list_commands = self.list_commands_for_help
        return super(SpecialHelpOrder, self).get_help(ctx)

    def list_commands_for_help(self, ctx):
        """reorder the list of commands when listing the help"""
        commands = super(SpecialHelpOrder, self).list_commands(ctx)
        return (c[1] for c in sorted(
            (self.help_priorities.get(command, 1), command)
            for command in commands))
    
    def command(self, *args, **kwargs):
        """Behaves the same as `click.Group.command()` except capture
        a priority for listing command names in help.
        """
        help_priority = kwargs.pop('help_priority', 1)
        help_priorities = self.help_priorities

        def decorator(f):
            cmd = super(SpecialHelpOrder, self).command(*args, **kwargs)(f)
            help_priorities[cmd.name] = help_priority
            return cmd

        return decorator


class CommandWrapper(object):

    def __init__(self, name, arguments=[], options={}, flags=[]):
        self.name = name
        self.arguments = arguments
        self.options = options
        self.flags = flags
        self.cmd_line = None

    def construct_command(self, option_prefix_char="-", flag_prefix_char="-"):
        _cmd_line = [self.name]
        for option, value in self.options.items():
            _cmd_line.extend([option_prefix_char + option, str(value)])
        for flag in self.flags:
            _cmd_line.append(flag_prefix_char + flag)
        self.cmd_line = _cmd_line
        return _cmd_line

    def run(self):
        if self.cmd_line:
            _logger.info("run {name} with command: {cmd}".format(name=self.name, cmd=" ".join(self.cmd_line)))
            #try:
            proc = subprocess.Popen(self.cmd_line, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            with proc.stdout as stdout:
                for line in iter(stdout.readline, b""):
                    _logger.debug(" [" + self.name +"]: " + line.decode('ASCII').strip())
            proc_exitcode = proc.wait()
            if proc_exitcode:
                _logger.error("run {name} failed with error code {err_code}\n".format(name=self.name, err_code=str(proc_exitcode)))
                raise subprocess.CalledProcessError(proc_exitcode, " ".join(self.cmd_line))


class CommandException(Exception):
    """Base class for exceptions in this module."""
    pass


def folder_exists(path_to_folder):
    """test path_to_folder is a folder and exists."""

    if os.path.isdir(path_to_folder) and os.path.exists(path_to_folder):
        return True
    else:
        return False


def create_directory(path_to_folder):
    """create folder when it doesn't exist."""

    if not folder_exists(path_to_folder):
        try:
            os.makedirs(path_to_folder)
        except Exception as e:
            if e.errno != errno.EEXIST:
                _logger.error(e)
            raise CommandException(e)


def get_prefix(input_file):
    """ make a file prefix based on input file
    :param input_file:
    :return: prefix of input_file
    """
    basename = os.path.basename(os.path.normpath(input_file))
    if "." in basename:
        prefix, _ = os.path.splitext(basename)
    else:
        prefix = basename
    _logger.info("no prefix is provided, '{0}' will be used as the output prefix".format(prefix))
    return prefix


def normalizer(input_freq_file, output_dir, prefix, scale_func=np.cbrt):
    """ 1) scale the data with numpy function
        2) normalize the data using Normalizer
    """

    output_norm_freq_file = os.path.join(output_dir, prefix+"_norm.tsv")


    freq_df = pd.read_csv(input_freq_file, sep="\t", index_col=0, header=0)
    freq_df.sort_index(inplace=True)

    # re-scale the freq df
    transform_freq_df = freq_df.apply(func=scale_func)

    # scale each value using sklearn.preprocessing.Normalizer()
    scaler = preprocessing.Normalizer()
    scaled_freq_df = scaler.fit_transform(transform_freq_df)
    scaled_freq_df = pd.DataFrame(scaled_freq_df, columns=transform_freq_df.columns, index=transform_freq_df.index)
    scaled_freq_df.to_csv(output_norm_freq_file, sep="\t", header=True, index=True, float_format='%.6f')

    return output_norm_freq_file
