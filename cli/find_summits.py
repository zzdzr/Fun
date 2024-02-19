import os
import cooler
import click
import logging

import numpy as np
import pandas as pd

from cli import cli
from lib.find_peaks import find_peak_prominence

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@cli.command()
@click.argument(
    "cool_path", metavar='COOL_PATH',
    type=str, nargs=1
)
@click.option(
    "--track",
    help = "The absolute path of SoN tracks for all chroms",
    type = str
)
@click.option(
    "--out_dir",
    help = "The absolute path for output directory",
    type = str
)

def generate_summits(cool_path, track, out_dir):
    """
    Find summits based on SoN

    """

    # load cooler and SoN track
    logger.info('Starting Summits detection...')
    clr = cooler.Cooler(cool_path)
    track_data = _load_track_data(track)

    # suffix
    resolution = clr.binsize
    suffix = f'_{resolution // 1000}kb.bed'

    out_dir = os.path.join(out_dir, 'SoN_summits/')
    os.makedirs(out_dir, exist_ok=True)

    for chrom in clr.chromnames:
        logger.info(f'Processing chromosome {chrom}...')
        _process_chromosome_data('chr' + chrom, track_data, out_dir, suffix)

    _merge_summits(out_dir, resolution)


def _load_track_data(track_path):
    """
    Load track data from the given path.

    """
    track_data = pd.read_table(
        track_path, header=None, sep='\t', names = ['chrom', 'start', 'end', 'SoN']
    )

    return track_data


def _process_chromosome_data(chrom, track, out_dir, suffix):
    """
    Process track data for a specific chromosome and write summits to file.

    """

    SoN_df = track[track.loc[:, 'chrom'] == chrom].copy()
    SoN_df.loc[:, 'SoN'] = SoN_df.loc[:, 'SoN'].clip(lower=0)

    # Get positions of summits
    poss, _ = find_peak_prominence(SoN_df.loc[:, 'SoN'].values)

    chr_cord_start = SoN_df.loc[:, 'start'].values[poss]
    chr_cord_end = SoN_df.loc[:, 'end'].values[poss]
    names = np.repeat('.', len(poss))
    strands = np.repeat('.', len(poss))

    # Get SoN scores
    SoN_values = SoN_df.loc[:, 'SoN'].values[poss]

    data_dict = {
        'chr': np.repeat(chrom, len(poss)),
        'start': chr_cord_start,
        'end': chr_cord_end,
        'name': names,
        'SoN': SoN_values,
        'strand': strands
    }

    df = pd.DataFrame(data_dict)
    df.to_csv(
        os.path.join(out_dir, chrom + suffix),
        sep='\t', header=True, index=None
    )

def _merge_summits(out_dir, resolution):
    """
    Merge individual summit files into one.

    """
    merged_file_name = f'Summits_{resolution}_merged.bed'
    merged_file_path = os.path.join(out_dir, merged_file_name)

    files = [f for f in os.listdir(out_dir) if f.endswith('.bed')]
    first = True

    if not os.path.exists(merged_file_path):
        with open(merged_file_path, 'w') as wfd:
            for f in files:
                file_path = os.path.join(out_dir, f)

                with open(file_path, 'r') as fd:
                    if first:
                        columns = next(fd)
                        wfd.write(columns)
                        first = False
                    else:
                        next(fd)

                    for line in fd:
                        wfd.write(line)

    else:
        logger.info('merged summits already exist, skip...')

