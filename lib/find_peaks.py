def find_peak_prominence(arr, max_dist=None):
    """Find the local maxima of an array and their prominence.
    The prominence of a peak is defined as the maximal difference between the
    height of the peak and the lowest point in the range until a higher peak.
    Parameters
    ----------
    arr : array_like
    max_dist : int
        If specified, the distance to the adjacent higher peaks is limited
        by `max_dist`.
    Returns
    -------
    loc_max_poss : numpy.array
        The positions of local maxima of a given array.
    proms : numpy.array
        The prominence of the detected maxima.
    """
    import numpy as np
    import warnings

    arr = np.asarray(arr)
    n = len(arr)
    max_dist = len(arr) if max_dist is None else int(max_dist)

    # Finding all local minima and maxima (i.e. points the are lower/higher than
    # both immediate non-nan neighbors).
    arr_nonans = arr[~np.isnan(arr)]
    idxs_nonans2idx = np.arange(arr.size)[~np.isnan(arr)]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)

        is_min_left = np.r_[False, arr_nonans[:-1] > arr_nonans[1:]]
        is_min_right = np.r_[arr_nonans[:-1] < arr_nonans[1:], False]
        is_loc_min = is_min_left & is_min_right
        loc_min_poss = np.where(is_loc_min)[0]
        loc_min_poss = idxs_nonans2idx[loc_min_poss]

        is_max_left = np.r_[False, arr_nonans[:-1] < arr_nonans[1:]]
        is_max_right = np.r_[arr_nonans[:-1] > arr_nonans[1:], False]
        is_loc_max = is_max_left & is_max_right
        loc_max_poss = np.where(is_loc_max)[0]
        loc_max_poss = idxs_nonans2idx[loc_max_poss]

    # For each maximum, find the position of a higher peak on the left and
    # on the right. If there are no higher peaks within the `max_dist` range,
    # just use the position `max_dist` away.
    left_maxs = -1 * np.ones(len(loc_max_poss), dtype=np.int64)
    right_maxs = -1 * np.ones(len(loc_max_poss), dtype=np.int64)

    for i, pos in enumerate(loc_max_poss):
        for j in range(pos - 1, -1, -1):
            if (arr[j] > arr[pos]) or (pos - j > max_dist):
                left_maxs[i] = j
                break

        for j in range(pos + 1, n):
            if (arr[j] > arr[pos]) or (j - pos > max_dist):
                right_maxs[i] = j
                break

    # Find the prominence of each peak with respect of the lowest point
    # between the peak and the adjacent higher peaks, on the left and the right
    # separately.
    left_max_proms = np.array(
        [
            (
                arr[pos] - np.nanmin(arr[left_maxs[i] : pos])
                if (left_maxs[i] >= 0)
                else np.nan
            )
            for i, pos in enumerate(loc_max_poss)
        ]
    )

    right_max_proms = np.array(
        [
            (
                arr[pos] - np.nanmin(arr[pos : right_maxs[i]])
                if (right_maxs[i] >= 0)
                else np.nan
            )
            for i, pos in enumerate(loc_max_poss)
        ]
    )

    # In 1D, the topographic definition of the prominence of a peak reduces to
    # the minimum of the left-side and right-side prominence.

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)

        max_proms = np.nanmin(np.vstack([left_max_proms, right_max_proms]), axis=0)

    # The global maximum, by definition, does not have higher peaks around it and
    # thus its prominence is explicitly defined with respect to the lowest local
    # minimum. This issue arises only if max_dist was not specified, otherwise
    # the prominence of the global maximum is already calculated with respect
    # to the lowest point within the `max_dist` range.
    # If no local minima are within the `max_dist` range, just use the
    # lowest point.
    global_max_mask = (left_maxs == -1) & (right_maxs == -1)
    if (global_max_mask).sum() > 0:
        global_max_idx = np.where(global_max_mask)[0][0]
        global_max_pos = loc_max_poss[global_max_idx]
        neighbor_loc_mins = (loc_min_poss >= global_max_pos - max_dist) & (
            loc_min_poss < global_max_pos + max_dist
        )
        if np.any(neighbor_loc_mins):
            max_proms[global_max_idx] = arr[global_max_pos] - np.nanmin(
                arr[loc_min_poss[neighbor_loc_mins]]
            )
        else:
            max_proms[global_max_idx] = arr[global_max_pos] - np.nanmin(
                arr[max(global_max_pos - max_dist, 0) : global_max_pos + max_dist]
            )

    return loc_max_poss, max_proms