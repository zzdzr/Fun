import bioframe
import cooler
import click
import pandas as pd
from lib.fountain_extension import *
from lib.generate_summits import *
from lib.merge_fountains import *
from lib.quality_filter import *
from lib.trans_to_bedpe import *
from lib.ks_test import *
from cli import cli

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@cli.command()
@click.argument(
    "cool_path", metavar='COOL_PATH',
    type=str, nargs=1
)
@click.option(
    "--half_width",
    help = "This parameter is related to the width of sampling and can be understood as follows: "
    "If the bin size is 10kb, then setting it to 2 means padding two 10kb bins on each side, "
    "thus extending the total to 50kb.",
    type = int
)
@click.option(
    "--ext_length",
    help = "This is actually the length of sampling box (bp). "
    "We recommand that 500Kb for 10kb resolution"
    "and 5Mb for 100kb resolution",
    type = int
)
@click.option(
    "--norm",
    help = "The normalization method for hic matrix"
    "(VC_SQRT, VC and KR normalization... based on your Cooler object)",
    default = 'VC_SQRT',
    show_default = True,
    type = str
)
@click.option(
    "--region_path",
    help = "Absolute path for fountains summits",
    type = str
)
@click.option(
    "--extension_pixels",
    help = "Array of locations we used to calculate dominance in extension. "
    "briefly, we just need to calculate the dominance at the specific pixel's position (for faster speed) "
    "This param needs further explanation in future version...",
    nargs = 3,
    type = int
)
@click.option(
    "--offset",
    help = "We do not consider length of extension below this threshold "
    "The term 'offset' refers to the distance perpendicular to the diagonal; "
    "pixels that are below this value in terms of their perpendicular distance "
    "from the diagonal do not participate in the computation (bp). ",
    default = 50000,
    show_default = True,
    type = int
)
@click.option(
    "--interval_length",
    help = "This param determines the length (bp) of sliding sheet containing "
    "multiple layers.",
    type = int
)
@click.option(
    "--coverage_ratio",
    help = "Given the targeted init bin for plumb,"
    "if the coverage of bins your collect is lower than "
    "this threshold, we consider it as NaNs."
    "For sparse matrix, you should consider this argument carefully",
    type = float
)
@click.option(
    "--output",
    help = "The absolute file path of output results",
    type = str
)
@click.option(
    "--p_value",
    help = "The threshold of p-value for K-S test",
    default=0.05,
    show_default=True,
    type = float
)
@click.option(
    "--signal_noise_background",
    help = "The threshold of SoN(fold change) for fountains",
    nargs=5,
    default=1.0,
    show_default=True,
    type = float
)
@click.option(
    "--max_merge_distance",
    help = "The maximum length we use to merge two close fountains",
    default=50000,
    show_default=True,
    type = int
)

def find_fountains(
    cool_path, half_width, ext_length,
    region_path, extension_pixels, offset,
    interval_length, coverage_ratio, output, norm=False,
    p_value = 0.05, signal_noise_background = 1.0,
    max_merge_distance = 20000,
):
    """
    Find fountains based on identified summits.

    """

    logger.info('Starting finding fountains...')

    # load Cooler object
    clr = cooler.Cooler(cool_path)

    # Determine resolution and load regions of summits
    resolution = clr.binsize
    region = bioframe.read_table(
        region_path, schema='bed'
    )

    # Bin array based on extension pixels (briefly, we just need to calculate the dominance at the specific pixel's position (for faster speed)）
    bin_array = np.arange(
        extension_pixels[0], extension_pixels[1], extension_pixels[2]
    )

    # Perform plumb calculation, get length and fold change of fountains
    logger.info('Calculate length of fountains and SoN (fold change)...')
    df = plumb(
        clr, half_width = half_width, extension_length=ext_length,
        norm=norm, regions=region, bin_array=bin_array, offset=offset,
        interval_length=interval_length, coverage_ratio = coverage_ratio
    )

    df = calculate_fountain_SoN(
        clr, df, norm=norm,
        half_width = half_width,
        coverage_ratio=coverage_ratio,
        offset = offset
    )

    # Perform K-S test
    logger.info('Perform K-S test...')
    df = make_ks_test(
        clr, regions=df, half_width = half_width,
        coverage_ratio=coverage_ratio, norm=norm,
        offset = offset
    )

    # Filter based on extension, p-value, and signal noise
    logger.info('Perform filter...')
    ext_bool = filter_extension(df)
    pvalue_bool = df['p_value'] < p_value

    for val in signal_noise_background:

        signal_bool = df['signal_noise_average_background'] > val
        df_tmp = df[(ext_bool & pvalue_bool & signal_bool)]
        df_tmp = merge_fountains(df_tmp.reset_index(drop=True), max_distance=max_merge_distance)
        output_tmp1 = f"{output}_{val}.tab"
        df_tmp.to_csv(output_tmp1, header=True, index=None, sep='\t')

        output_tmp2 = f"{output}_{val}.bedpe"
        bedpe = dataframe_to_bedpe(df_tmp, resolution=resolution)
        bedpe.to_csv(output_tmp2, header=True, index=None, sep='\t')

    logger.info('Complete!')
