import copy
import numba
import cooler
import warnings
import numpy as np
import matplotlib.pyplot as plt

@numba.jit(nopython=True)
def plumb_numbda(half_width, init_bin, extension_length, n_bins, offset, resolution):
    '''

    Notes:
    We use numba to accelerate the loops, therefore,
    we could fastly generate the coordinates we used to pileup

    You should mention that the init_bin should not exceed the boundary of matrix,
    it will result some unknown errors.

    * If the extension length < resolution, we only obtain the initial value for the first pixel (init pixel)
    * If exceed the boundary, the element is empty
    '''

    if np.isnan(extension_length):
        raise ValueError('nan extension length')

    extension_bins = (extension_length // resolution) * 2

    if offset != 0:
        offset_bins = (offset // resolution) * 2 + 1
    else:
        offset_bins = 0

    init_row_idx, init_col_idx = np.arange(init_bin - half_width, init_bin + half_width + 1), \
                                 np.arange(init_bin - half_width, init_bin + half_width + 1)

    assert offset_bins < extension_bins+1, \
        'The offset should not exceed length of extension'

    for i in range(offset_bins, extension_bins + 1):
        if i % 2 == 0:
            row_idx = init_row_idx - i // 2
            col_idx = init_col_idx + i // 2

        else:
            row_idx = row_idx[:-1] if len(row_idx) > 1 else row_idx
            col_idx = col_idx[:-1] + 1 if len(col_idx) > 1 else col_idx + 1

        # Do not exceed the boundaries
        cond = (row_idx >= 0) & (col_idx <= n_bins - 1)
        row_idx = row_idx[cond]
        col_idx = col_idx[cond]

        yield row_idx, col_idx


@numba.jit(nopython=True)
def nb_sum(x):
    sum_ = 0
    for i in range(len(x)):
        sum_ += x[i]
    return sum_

def plumb_sum(mat, half_width, init_bin, extension_length, resolution, offset):
    '''
    calculate sum of interactions in sampling box at target bin

    Parameters
    ----------
    mat: ndarray object
        hic matrix

    half_width: int object
        number of bins of padding regions for its center

    extension_length: int object
        the length of extension (bp)

    Returns
    -------
    val : float object
        sum of interaction at each distance
    '''

    for row_idx, col_idx in plumb_numbda(
        half_width=half_width,
        init_bin=init_bin,
        extension_length=extension_length,
        n_bins=mat.shape[0],
        offset=offset,
        resolution=resolution
    ):

        if row_idx.size == 0 or col_idx.size == 0:
            val = np.nan
        else:
            val = np.nansum(mat[row_idx, col_idx])

        yield val

def global_diag_plumb_sum(mat, half_width, extension_length, offset):
    '''
    This function can perform "plumb_sum" at all bins along the matrix

    '''
    for _ in range(n_bins):
        diag_plumb_val = np.asarray([
            i for i in plumb_sum(mat, half_width, _, extension_length, offset)
        ])

        yield diag_plumb_val

def plumb_layer_values(mat, half_width, init_bin, extension_length, resolution, offset):
    '''
    calculate sum of interactions in sampling box at target bin

    Parameters
    ----------
    mat: ndarray object
        hic matrix

    half_width: int object
        number of bins of padding regions for its center

    extension_length: int object
        the length of extension (bp)

    Returns
    -------
    val : float object
        sum of interaction at each distance
    '''

    for row_idx, col_idx in plumb_numbda(
        half_width=half_width,
        init_bin=init_bin,
        extension_length=extension_length,
        n_bins=mat.shape[0],
        offset=offset,
        resolution=resolution
    ):

        layers = []
        if row_idx.size == 0 or col_idx.size == 0:
            val = np.nan
        else:
            val = mat[row_idx, col_idx]

        layers.append(val)

    return layers


def plumb_mean(mat, half_width, init_bin, extension_length, resolution, offset):
    '''
    calculate sum of interactions in sampling box at target bin

    Parameters
    ----------
    mat: ndarray object
        hic matrix

    half_width: int object
        number of bins of padding regions for its center

    extension_length: int object
        the length of extension (bp)

    Returns
    -------
    val : float object
        sum of interaction at each distance
    '''

    for row_idx, col_idx in plumb_numbda(
        half_width=half_width,
        init_bin=init_bin,
        extension_length=extension_length,
        n_bins=mat.shape[0],
        offset=offset,
        resolution=resolution
    ):

        if row_idx.size == 0 or col_idx.size == 0:
            val = np.nan
        else:
            val = np.nanmean(mat[row_idx, col_idx])

        yield val