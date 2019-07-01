# -*- coding: utf-8 -*-


import os
import sys
import click
import logging


_logger = logging.getLogger(__name__)


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
    
    def add_command(self, *args, **kwargs):
        """overwrite the add_command method to allow help_priority parameter"""
        return self.command(*args, **kwargs)