import bioframe
from cli import cli
from lib.signal_over_noise import *
from lib.util import *

import click
import os
import logging

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
    default = 0.2,
    show_default = True,
    type = float
)
@click.option(
    "--ext_length",
    help = "Extension length for given init bin"
    "We recommand that 500Kb for 10kb resolution"
    "and 5Mb for 100kb resolution",
    default = 500000,
    show_default = True,
    type = int
)
@click.option(
    "--padding_width",
    help = "N of bins for given init bin. "
    "For a given sampling box, the width is (2N+1)*resolution",
    default = 2,
    show_default = True,
    type = int
)
@click.option(
    "--offset",
    help = "We do not consider length of extension below this threshold",
    default = 20000,
    show_default = True,
    type = int
)
@click.option(
    "--merge",
    help = "Whether merge .bedgraph track for all chromosmes",
    default = True,
    show_default = True,
    type = bool
)
@click.option(
    "--use_mean",
    help = "Use mean else median",
    default = False,
    show_default = True,
    type = bool
)
@click.option(
    "--chromsize_path",
    help = "file containing chromsize",
    type = str
)
@click.option(
    "--bedGraphToBigWig",
    help = "path to your bedGraphToBigWig from UCSC",
    type = str
)

def calculate_SoN_score(
    cool_path, out_dir, chromsize_path, norm = False, coverage_ratio = 0.2,
    ext_length = 500000, padding_width = 2, offset=20000, merge = True, use_mean = False
):
    """
    Calculate signal-over-noise (SoN) score for all chromosomes.

    Args:
        cool_path (str): Path to the cooler file.
        out_dir (str): Output directory for SoN tracks.
        chromsize_path (str): Path to the chromosome size file.
        norm (bool): Normalization strategy.
        coverage_ratio (float): Coverage ratio.
        ext_length (int): Extension length.
        padding_width (int): Half-width of padding.
        offset (int): Offset value.
        merge (bool): Whether to merge the generated bedgraph files.
        use_mean (bool): Use mean for SoN score calculation.

    """
    logger.info('Starting SoN calculation...')


    # Load cooler object
    logger.info(f'Loading cooler object from {cool_path}...')
    clr = cooler.Cooler(cool_path)
    resolution = clr.binsize

    # load chromsize
    logger.info(f'Loading chromosome sizes from {chromsize_path}...')
    chromsize = bioframe.read_chromsizes(
        chromsize_path, natsort=True
    )

    # Make directory for SoN tracks
    out_dir = os.path.join(out_dir, f'SoN_track_{resolution}')
    os.makedirs(out_dir, exist_ok=True)
    logger.info(f'Created directory for SoN tracks: {out_dir}')

    logger.info(f'Your normalization strategy is {norm}')
    for chrom in clr.chromnames:
        logger.info(f'Processing chromosome {chrom}...')

        # Output path
        output = os.path.join(
            out_dir, f'chr{chrom}_{resolution}_SoN.bedgraph'
        )

        # Ignore if file already exists
        if os.path.exists(output):
            continue

        mat = clr.matrix(balance=norm).fetch(chrom)

        # Calculate SoN score
        SoN_score = calculate_signal_noise_ratio_score(
            mat=mat, extension_length=ext_length,
            resolution=resolution, half_width= padding_width,
            offset=offset, coverage_ratio=coverage_ratio, use_mean = use_mean
        )

        # NaN or infinite values are regarded as 0
        SoN_score = np.asarray(list(SoN_score))
        SoN_score[np.isnan(SoN_score)] = 0
        SoN_score[np.isinf(SoN_score)] = 0

        # Get coordinates of SoN score
        chr_cord_start = [i * resolution for i in range(len(SoN_score))]
        chr_cord_end = [(i + 1) * resolution for i in range(len(SoN_score))]

        data_dict = {
            'chr': np.repeat('chr' + chrom, len(chr_cord_start)),
            'start': chr_cord_start,
            'end': chr_cord_end,
            'value': SoN_score
        }


        SoN_track = align_track_with_chromsize(
            pd.DataFrame(data_dict), chromsize
        )

        SoN_track.to_csv(
            output, sep='\t', header=None,index=None
        )

    if merge:
        merged_file_path = os.path.join(
            out_dir, f'SoN_{resolution}_merged.bedgraph'

        )
        os.system(
            'cat %s/*.bedgraph | sort -k1,1 -k2,2n > %s' % (out_dir, merged_file_path)
        )

        bigwig_output_path = os.path.join(
            out_dir, f'SoN_{resolution}_merged.bw'
        )

        bedgraph_to_bigwig(
            bedGraphToBigWig, merged_file_path, chromsize_path, bigwig_output_path
        )

    logger.info('SoN score calculation completed.')