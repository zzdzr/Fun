# -*- coding: utf-8 -*-
import click
import sys

CONTEXT_SETTING = {
    'help_option_names':['-h', '--help']
}

@click.version_option('1.0', '-V', '--version')
@click.group(context_settings = CONTEXT_SETTING)
@click.option('-v', '--verbose', help = 'Verbose test', is_flag=True, default=False)
def cli(verbose):
    if verbose:
        print('Is verbose')

from . import (
    calculate_extension_infor,
    calculate_SoN,
    find_summits
)

