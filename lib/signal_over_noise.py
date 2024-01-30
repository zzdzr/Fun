from multiprocessing import Pool
from scipy.stats import ks_2samp
from .diagonal_plumb import *
import pandas as pd
import numpy as np
import cooler
import sys
import os
import scipy


def ks_test(mat, half_width, mid, extension_length, offset, coverage_ratio, resolution):

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


def calculate_coverage(plumb_signal):
    '''
    Calculate coverage ratio when perform plumb

    '''

    plumb_signal[np.isnan(plumb_signal)] = 0

    # count ratio of non-zero element
    cov_ratio = np.count_nonzero(plumb_signal) / len(plumb_signal)

    return cov_ratio


def set_background(
    mat, half_width, extension_length, bin_idx, offset, resolution,
    orientation = None, coverage_ratio = None
):
    '''
    Used to calculate the background of plumbs, which will help
    to calculate the fountainess future

    Notes:
    The function is directly based on the plumb-pileup function

    '''
    if coverage_ratio is None:
        raise ValueError('Empty coverage ratio.')

    # Upstream background
    if orientation == 'upstream':
        diag_plumb_val = np.asarray(
            [i for i in plumb_sum(
                mat = mat, half_width = half_width,
                init_bin = bin_idx - 2*half_width-1,
                extension_length = extension_length,
                offset = offset, resolution = resolution
            )]
        )

        if len(diag_plumb_val) >= 1:
            cov_ratio = calculate_coverage(diag_plumb_val)

        else:
            cov_ratio = 0

        if cov_ratio <= coverage_ratio:
            diag_plumb_val = np.nan


    elif orientation == 'downstream':
        diag_plumb_val = np.asarray(
            [i for i in plumb_sum(
                mat = mat, half_width = half_width,
                init_bin = bin_idx + 2*half_width+1,
                extension_length = extension_length,
                offset = offset, resolution = resolution
            )]
        )
        if len(diag_plumb_val) >= 1:
            cov_ratio = calculate_coverage(diag_plumb_val)

        else:
            cov_ratio = 0

        if cov_ratio <= coverage_ratio:
            diag_plumb_val = np.nan

    elif orientation is None:
        raise TypeError('Please input orientation of background')

    else:
        raise TypeError('Invalid orientation')

    return diag_plumb_val


def center_area_plumb_sum(
    mat, half_width, extension_length,
    bin_idx, offset, resolution, coverage_ratio = None
):

    if coverage_ratio is None:
        raise ValueError('Empty coverage ratio.')

    target_plumb_val = np.asarray(
        [i for i in plumb_sum(
            mat = mat,
            half_width = half_width,
            init_bin = bin_idx, extension_length = extension_length,
            resolution = resolution, offset = offset
        )]
    )

    if len(target_plumb_val) >= 1:
        cov_ratio = calculate_coverage(target_plumb_val)

    else:
        cov_ratio = 0

    if cov_ratio <= coverage_ratio:
        target_plumb_val = np.nan


    return target_plumb_val


def calculate_signal_noise_ratio_score(
    mat, half_width, extension_length,
    resolution, offset, coverage_ratio = 0.2,
    use_mean = True
):
    '''
    calculate a statistics to measure the strength of fountain

    Parameters
    ----------
    mat: ndarray object
        hic matrix

    half_width: int object
        number of bins of padding regions for its center

    extension_length: int object
        Length of extension (bp)

    Returns
    -------
    fountain score: float object
        a kind of statistics to measure the strength of fountain

    '''

    for idx in range(mat.shape[0]):
        up_bkg = set_background(
            mat=mat, half_width=half_width, extension_length=extension_length,
            offset=offset, bin_idx=idx, orientation = 'upstream',
            coverage_ratio=coverage_ratio, resolution=resolution
        )

        down_bkg = set_background(
            mat=mat, half_width=half_width, extension_length=extension_length,
            offset=offset, bin_idx=idx, orientation = 'downstream',
            coverage_ratio=coverage_ratio, resolution=resolution
        )

        tkg_plumb = center_area_plumb_sum(
            mat=mat, half_width=half_width, extension_length=extension_length,
            offset=offset, bin_idx=idx, coverage_ratio=coverage_ratio, resolution=resolution
        )

        if np.all(
            [isinstance(i, np.ndarray) for i in [up_bkg, down_bkg, tkg_plumb]]
        ):
            # gradient for upstream / downstream
            bkg_ave = np.nanmean(
                [np.nansum(up_bkg), np.nansum(down_bkg)]
            )

            bkg_gradient = np.nansum(tkg_plumb) / bkg_ave

            # set np.nan/inf to 0
            if np.isnan(bkg_gradient) or np.isinf(bkg_gradient):
                bkg_gradient = 1

            if use_mean:
                fountain_score = np.nanmean(tkg_plumb) * np.log(bkg_gradient)
            else:
                fountain_score = np.nanmedian(tkg_plumb) * np.log(bkg_gradient)

        else:
            fountain_score = np.nan

        yield fountain_score


def worker(args):
    mat, half_width, extension_length, resolution, offset, coverage_ratio, idx = args
    tkg_plumb = center_area_plumb_sum(
        mat, half_width, extension_length, idx, offset, resolution, coverage_ratio
    )

    print(f"This is bin {idx}")

    if isinstance(tkg_plumb, np.ndarray):
        return np.nansum(tkg_plumb)
    else:
        return np.nan

def calculate_strength_parallel(mat, half_width, extension_length, resolution, offset, coverage_ratio, CPU):

    print('Calculating parallel......')
    print(f'------Your input CPU is {CPU}------')

    with Pool(CPU) as pool:
        args_list = [(mat, half_width, extension_length, resolution, offset, coverage_ratio, idx) for idx in range(mat.shape[0])]
        results = pool.map(worker, args_list)
    return results


def calculate_strength_single(mat, half_width, extension_length, resolution, offset, coverage_ratio):

    for idx in range(mat.shape[0]):

        tkg_plumb = center_area_plumb_sum(
            mat, half_width, extension_length, idx, offset, resolution, coverage_ratio
        )

        if isinstance(tkg_plumb, np.ndarray):
            yield np.nansum(tkg_plumb)
        else:
            yield np.nan



def read_chrom_sizes(chrominfo_file, prefix=""):
    with open(chrominfo_file, 'r') as f:
        chrom_sizes = {}
        for line in f:
            chrom, size = line.strip().split()
            chrom_id = f"{prefix}{chrom.replace('chr', '')}"
            chrom_sizes[chrom_id] = int(size)
    return chrom_sizes


