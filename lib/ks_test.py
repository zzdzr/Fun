# make offset = 0
import numpy as np
from Replihic.lib.fountain_extension import *

def ks_test(mat, half_width, mid, extension_length, coverage_ratio, resolution, offset):

    up_bkg_per_dist = set_background(
        mat=mat, half_width=half_width, extension_length=extension_length, resolution=resolution,
        offset=offset, bin_idx=mid, orientation='upstream', coverage_ratio=coverage_ratio
    )

    down_bkg_per_dist = set_background(
        mat=mat, half_width=half_width, extension_length=extension_length, resolution = resolution,
        offset=offset, bin_idx=mid, orientation='downstream', coverage_ratio=coverage_ratio
    )

    tkg_plumb_per_dist = center_area_plumb_sum(
        mat=mat, half_width=half_width, extension_length=extension_length, resolution=resolution,
        offset=offset, bin_idx=mid, coverage_ratio=coverage_ratio
    )

    if np.all(
        [isinstance(i, np.ndarray) for i in [up_bkg_per_dist, down_bkg_per_dist, tkg_plumb_per_dist]]
    ):
        bkg_ave = np.nanmean([up_bkg_per_dist, down_bkg_per_dist], axis=0)
        kstest_res = ks_2samp(tkg_plumb_per_dist, bkg_ave)

        return kstest_res

    else:
        return np.nan

# make offset = 0
def make_ks_test(clr, clr_control, regions, half_width, norm, coverage_ratio, offset, use_control=True):
    '''
    for each of candidate fountain, we use Kolmogorov–Smirnov test to verify its credibility

    Parameters
    ----------
    mat: ndarray object


    Returns
    -------
    regions: DataFrame object
    regions of candidates containing p-value and FDR

    '''
    regions = regions.copy()
    resolution = clr.binsize

    if not 'max_extension' in regions.columns:
        raise TypeError(
            'Invalid dataframe, you have to calculate the length of extension'
        )

    chrom_list = []
    ks_res_list = []
    for index, row in regions.iterrows():
        chrom = row['chrom'][3:]
        mid = (row['end'] + row['start']) // resolution // 2
        extension_length = row['max_extension'] * 1000

        if np.isnan(extension_length) or offset >= extension_length:
            ks_res = np.nan

        else:
            if chrom in chrom_list:
                if use_control:
                    ks_res = ks_test(
                        mat=ratio_mat, half_width=half_width, mid=mid, extension_length=extension_length,
                        coverage_ratio=coverage_ratio, resolution=resolution, offset=offset
                    )
                else:
                    ks_res = ks_test(
                        mat=mat, half_width=half_width, mid=mid, extension_length=extension_length,
                        coverage_ratio=coverage_ratio, resolution=resolution, offset=offset
                    )

            else:
                chrom_list.append(chrom)
                mat = clr.matrix(balance = norm).fetch(chrom)

                if use_control:
                    mat_control = clr_control.matrix(balance = norm).fetch(chrom)
                    ratio_mat = np.divide(mat, mat_control, out=np.zeros_like(mat), where=mat_control != 0)
                    ratio_mat = np.nan_to_num(ratio_mat, nan=0)
                    ks_res = ks_test(
                        mat=ratio_mat, half_width=half_width, mid=mid, extension_length=extension_length,
                        coverage_ratio=coverage_ratio, resolution=resolution, offset=offset
                    )
                else:
                    ks_res = ks_test(
                        mat=mat, half_width=half_width, mid=mid, extension_length=extension_length,
                        coverage_ratio=coverage_ratio, resolution=resolution, offset=offset
                    )

        ks_res_list.append(ks_res)

    regions['p_value'] = [
        i.pvalue if isinstance(
            i, scipy.stats._stats_py.KstestResult
        ) else np.nan for i in ks_res_list
    ]

    return regions
