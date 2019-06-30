# -*- coding: utf-8 -*-

"""Console script for bears."""
import sys
import click


# check here for custom click group and command:
# https://github.com/sphinx-contrib/sphinxcontrib-versioning/blob/master/sphinxcontrib/versioning/__main__.py

# change the order of the subcommands:
# https://stackoverflow.com/questions/47972638/how-can-i-define-the-order-of-click-sub-commands-in-help/47984810#47984810

@click.group()
@click.argument('assembly', type=str)
@click.argument('bam_files', nargs=-1)
@click.pass_context
def main(ctx, args=None):
    """Console script for bears."""
    click.echo("Replace this message by putting your code into "
               "bears.cli.main")
    click.echo("See click documentation at http://click.pocoo.org/")
    return 0


@main.command()
@click.pass_context
def profile(ctx):
    print('input_assembly is {assembly}'.format(assembly = ctx.obj['assembly']))


@main.command()
@click.pass_context
def bin(ctx):
    print('input_assembly is {assembly}'.format(assembly = ctx.obj['assembly']))


if __name__ == "__main__":
    sys.exit(main()) 
