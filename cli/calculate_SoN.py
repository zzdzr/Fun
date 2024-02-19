import bioframe
import logging
import glob
import click
import os

from cli import cli
from lib.signal_over_noise import *
from lib.util import *



logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@cli.command()
@click.argument(
    "cool_path", metavar = 'COOL_PATH',
    type = str, nargs = 1
)
@click.option(
    "--out_dir",
    help = "The output directory of SoNs",
    type = str
)
@click.option(
    "--norm",
    help = "The normalization method for hic matrix"
    "(VC_SQRT, VC and KR normalization)",
    default = 'VC_SQRT',
    show_default = True,
    type = str
)
@click.option(
    "--coverage_ratio",
    help = "Given the targeted init bin for plumb,"
    "if the coverage of bins your collect is lower than "
    "this threshold, we consider it as NaNs."
    "For sparse matrix, you should consider this argument carefully",
    default = 0,
    show_default = True,
    type = float
)
@click.option(
    "--ext_length",
    help = "This is actually the length of sampling box (bp). "
    "We recommand that 500Kb for 10kb resolution"
    "and 5Mb for 100kb resolution",
    default = 500000,
    show_default = True,
    type = int
)
@click.option(
    "--padding_width",
    help = "This parameter is related to the width of sampling and can be understood as follows: "
    "If the bin size is 10kb, then setting it to 2 means padding two 10kb bins on each side, "
    "thus extending the total to 50kb.",
    default = 2,
    show_default = True,
    type = int
)
@click.option(
    "--offset",
    help = "We do not consider length of extension below this threshold "
    "The term 'offset' refers to the distance perpendicular to the diagonal; "
    "pixels that are below this value in terms of their perpendicular distance "
    "from the diagonal do not participate in the computation (bp). ",
    default = 20000,
    show_default = True,
    type = int
)
@click.option(
    "--integrate",
    help = "Whether merge .bedgraph track for all chromosmes",
    default = True,
    show_default = True,
    type = bool
)
@click.option(
    "--use_mean",
    help = "Use mean for SoN score calculation in pixels from sampling regions else median",
    default = False,
    show_default = True,
    type = bool
)
@click.option(
    "--chromsize_path",
    help = "file containing chromsize",
    type = str
)


def calculate_SoN_score(
    cool_path, out_dir, chromsize_path, norm=False,
    coverage_ratio=0.2, ext_length=500000, padding_width=2,
    offset=20000, integrate=True, use_mean=False
):
    """
    Calculate signal-over-noise (SoN) score for a specific chromosome.

    Returns:
        Dataframe containing SoN scores and corresponding genomic coordinates.

    """
    logger.info('Starting SoN calculation...')
    clr, resolution = _load_cooler_object(cool_path)
    chromsize = bioframe.read_chromsizes(chromsize_path, natsort=True)
    out_dir_path = _create_output_directory(out_dir, resolution)

    for chrom in clr.chromnames:

        logger.info(f'Processing chromosome {chrom}...')

        # calculate SoN score
        SoN_track = _calculate_SoN_score(
            clr, chrom, norm, ext_length, padding_width,
            offset, coverage_ratio, use_mean, chromsize, out_dir_path, resolution
        )

        # output tracks
        _write_SoN_tracks(
            out_dir_path, chrom, resolution, SoN_track
        )

    logger.info('SoN score calculation completed.')

    if integrate:
        logger.info('Merging all chromosome tracks into one file...')
        _merge_bedgraph_files(out_dir_path, resolution)


def _load_cooler_object(cool_path):
    """
    Load cooler object from the given path.

    Args:
        cool_path (str): Path to the cooler file.

    Returns:
        Cooler object and its resolution.
    """
    logger.info(f'Loading cooler object from {cool_path}...')
    clr = cooler.Cooler(cool_path)
    resolution = clr.binsize
    return clr, resolution


def _create_output_directory(out_dir, resolution):
    """
    Create output directory for SoN tracks.

    Args:
        out_dir (str): Base output directory.
        resolution (int): Resolution of the cooler file.

    Returns:
        Path to the created output directory.
    """
    out_dir_path = os.path.join(out_dir, f'SoN_track_{resolution}')
    os.makedirs(out_dir_path, exist_ok=True)
    logger.info(f'Created directory for SoN tracks: {out_dir_path}')
    return out_dir_path


def _calculate_SoN_score(clr, chrom, norm, ext_length, padding_width, offset, coverage_ratio, use_mean, chromsize, out_dir_path, resolution):
    """
    Calculate signal-over-noise (SoN) score for a specific chromosome.

    """

    # Check if SoN track already exists
    output_path = os.path.join(out_dir_path, f'chr{chrom}_{resolution}_SoN.bedgraph')
    if os.path.exists(output_path):
        logger.info(f'SoN track for chromosome {chrom} already exists. Skipping...')

        return None

    mat = clr.matrix(balance=norm).fetch(chrom)
    SoN_score = calculate_signal_noise_ratio_score(
        mat=mat, extension_length=ext_length,
        resolution=resolution, half_width=padding_width,
        offset=offset, coverage_ratio=coverage_ratio, use_mean=use_mean
    )

    # Handle NaN or infinite values
    SoN_score = np.asarray(list(SoN_score))
    SoN_score = np.nan_to_num(SoN_score, nan = 0, posinf=0, neginf=0)

    # Get genomic coordinates
    chr_cord_start = [i * clr.binsize for i in range(len(SoN_score))]
    chr_cord_end = [(i + 1) * clr.binsize for i in range(len(SoN_score))]
    data_dict = {
        'chr': np.repeat('chr' + chrom, len(chr_cord_start)),
        'start': chr_cord_start,
        'end': chr_cord_end,
        'value': SoN_score
    }

    SoN_track = align_track_with_chromsize(pd.DataFrame(data_dict), chromsize)
    return SoN_track

def _write_SoN_tracks(out_dir, chrom, resolution, SoN_track):
    """
    Write SoN tracks to a bedgraph file.

    Args:
        out_dir (str): Output directory for SoN tracks.
        chrom (str): Chromosome name.
        resolution (int): Resolution of the cooler file.
        SoN_track (DataFrame): Dataframe containing SoN scores and genomic coordinates.
    """
    output_path = os.path.join(out_dir, f'chr{chrom}_{resolution}_SoN.bedgraph')

    if SoN_track is not None:

        SoN_track.to_csv(
            output_path, sep='\t', header=None, index=None
        )

        logger.info(f'SoN track for chromosome {chrom} written to {output_path}')

def _merge_bedgraph_files(out_dir, resolution):
    """
    Merge all .bedgraph files into one.

    Args:
        out_dir (str): The directory where the .bedgraph files are stored.
        resolution (int): Resolution of the cooler file.
    """
    merged_file_name = f'SoN_{resolution}_merged.bedgraph'
    merged_file_path = os.path.join(out_dir, merged_file_name)

    # Find all bedgraph files in the output directory
    bedgraph_files = glob.glob(os.path.join(out_dir, '*.bedgraph'))

    # Read all bedgraph files and concatenate them
    df_list = [
        pd.read_csv(file, sep='\t', header=None) for file in bedgraph_files
    ]

    merged_df = pd.concat(df_list, axis=0)

    # Sort by chromosome and start position
    merged_df.sort_values(by=[0, 1], inplace=True)

    # Write merged data to a new bedgraph file
    merged_df.to_csv(merged_file_path, sep='\t', header=False, index=False)

    logger.info(f'Merged bedgraph file created at {merged_file_path}')