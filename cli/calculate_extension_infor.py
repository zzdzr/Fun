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


@cli.command()
@click.argument(
    "cool_path", metavar='COOL_PATH',
    type=str, nargs=1
)
@click.option(
    "--half_width",
    help = "The number of bins that pad summit "
    "This param can decide the final width of sampling box((2*half_width+1) * resolution) ",
    type = int
)
@click.option(
    "--ext_length",
    help = "Extension length for given init bin. "
    "We recommand that 500Kb for 10kb resolution "
    "and 5Mb for 100kb resolution",
    type = int
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
    "--region_path",
    help = "All fountain summits regions",
    type = str
)
@click.option(
    "--extension_pixels",
    help = "Array of locations we used to calculate dominance in extension",
    nargs = 3,
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
    "--interval_length",
    help = "The length of sampling box in extension",
    type = int
)
@click.option(
    "--coverage_ratio",
    help = "When plumb/extension, if the number of non-nan bins / total < threshold "
    "for this summit, the plumb values will be nan.",
    type = float
)
@click.option(
    "--output",
    help = "The absolute path with prefix of output result",
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
    help = "The threshold of SoN for fountain",
    nargs=5,
    default=1.0,
    show_default=True,
    type = float
)
@click.option(
    "--max_merge_distance",
    help = "The maximum length we use to merge two close fountains",
    default=20000,
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

    # load cooler object
    clr = cooler.Cooler(cool_path)

    resolution = clr.binsize
    region = bioframe.read_table(
        region_path, schema='bed'
    )

    # get bin array
    bin_array = np.arange(
        extension_pixels[0], extension_pixels[1], extension_pixels[2]
    )

    # calculate maximum extension length
    df = plumb(
        clr, half_width = half_width, extension_length=ext_length,
        norm=norm, regions=region, bin_array=bin_array, offset=offset,
        interval_length=interval_length, coverage_ratio = coverage_ratio
    )

    # calculate median strength of fountain and its SoN score
    df = calculate_fountain_SoN(
        clr, df, norm=norm,
        half_width = half_width,
        coverage_ratio=coverage_ratio,
        offset = offset
    )

    df = make_ks_test(
        clr, regions=df, half_width = half_width,
        coverage_ratio=coverage_ratio, norm=norm,
        offset = offset
    )

    ext_bool = filter_extension(df)
    pvalue_bool = df['p_value'] < p_value

    for val in signal_noise_background:

        signal_bool = df['signal_noise_average_background'] > val

        df_tmp = df[(ext_bool & pvalue_bool & signal_bool)]
        df_tmp = merge_fountains(df_tmp.reset_index(drop=True), max_distance=max_merge_distance)
        output_tmp1 = output + '_{}.bed'.format(val)
        df_tmp.to_csv(output_tmp1, header=True, index=None, sep='\t')

        output_tmp2 = output + '_{}.bedpe'.format(val)
        bedpe = dataframe_to_bedpe(df_tmp, resolution=resolution)
        bedpe.to_csv(output_tmp2, header=True, index=None, sep='\t')
