from .signal_over_noise import *
import logging

# logging.basicConfig(format='%(levelname)s:%(funcName)s:%(message)s', level=logging.DEBUG)
def find_maximum_extension(signal_track, resolution, interval_length = 50000, threshold = 0.5):

    if not isinstance(signal_track, np.ndarray):
        raise TypeError('Invalid signal track, please use ndarray!')

    n_bins = (interval_length // resolution) * 2 - 1


    #half_width = interval_length // 2 // resolution * 2
    half_width = interval_length//resolution - 1

    for i in range(0, len(signal_track)):

        if i < half_width or i + half_width > len(signal_track) - 1:
            n_dominance = np.nan

        else:
            signal = signal_track[i - half_width: i + half_width + 1]
            n_dominance = signal[signal <= 0].size

        # we could stop at the bin i
        if n_dominance / n_bins >= threshold:

            # then we should find the bin with max value
            #max_idx = np.argmax(signal_track[i - half_width: i+1])

            max_idx = np.argmax(signal_track[i-half_width:i+half_width+1])

            obtained_idx = i - n_bins // 2 + max_idx

            return obtained_idx

    return n_bins


def calculate_dominance(signal_track, bin_array, threshold = 0):

    if bin_array[0] < 1:
        raise ValueError('The first index should greater than 0, otherwise empty')

    if isinstance(signal_track, np.ndarray):
        percent_list = []
        # logging.debug('The shape of signal_track is %d' % signal_track.shape)
        # logging.debug('The length of bin array is %d' % bin_array.shape)

        if bin_array.size > signal_track.size:
            raise ValueError(
                'The index in bin_array should not exceed the bound of ndarray'
            )

        for idx in bin_array:
            signal = signal_track[:idx]
            percent = signal[signal > threshold].size / signal.size
            percent_list.append(percent)

    else:
        raise TypeError('Signal track must be type of ndarray!')

    percent_list = np.asarray(percent_list)

    if bin_array[0] > 1:
        percent_list = np.append(np.repeat(np.nan, bin_array[0] - 1), percent_list)

    return percent_list

def center_background_plumb(
    mat, half_width, extension_length, init_bin, resolution,
    offset, bin_array, coverage_ratio = 0.2, interval_length = 50000, threshold = 0.5
):

    tkg_plumb, up_bkg, down_bkg = set_sampling_box(
        mat = mat, half_width = half_width, extension_length = extension_length,
        offset = offset, init_bin = init_bin, coverage_ratio=coverage_ratio, resolution=resolution
    )

    if np.all(
        [isinstance(i, np.ndarray) for i in [up_bkg, down_bkg, tkg_plumb]]
    ):

        bkg_ave = np.nanmean(
            [up_bkg, down_bkg], axis = 0
        )

        tkg_bkg_subtract = tkg_plumb - bkg_ave
        maximum_extension = find_maximum_extension(
            tkg_bkg_subtract, resolution = resolution,
            interval_length = interval_length, threshold = threshold
        )

        percent_list = calculate_dominance(
             tkg_bkg_subtract, bin_array = bin_array
        )

    else:
        return np.nan, np.nan

    return maximum_extension, percent_list

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

def set_sampling_box(
    mat, half_width, extension_length, init_bin, resolution,
    offset, coverage_ratio = 0.2
):

    tkg_plumb = center_area_plumb_sum(
        mat=mat, half_width=half_width, extension_length=extension_length,
        offset=offset, bin_idx=init_bin, coverage_ratio=coverage_ratio, resolution=resolution
    )

    up_bkg = set_background(
        mat=mat, half_width=half_width, extension_length=extension_length,
        offset=offset, bin_idx=init_bin, orientation='upstream',
        coverage_ratio=coverage_ratio, resolution=resolution
    )

    down_bkg = set_background(
        mat=mat, half_width=half_width, extension_length=extension_length,
        offset=offset, bin_idx=init_bin, orientation='downstream',
        coverage_ratio=coverage_ratio, resolution=resolution
    )

    return tkg_plumb, up_bkg, down_bkg

def calculate_fountain_SoN(
    clr, regions, half_width,
    coverage_ratio, norm, offset
):

    regions = regions.copy()
    regions = regions.sort_values(by = 'chrom').reset_index(drop = True)

    if not set(
        ['chrom', 'start', 'end', 'max_extension']
    ).issubset(regions.columns):

        raise TypeError(
            'Invalid dataframe \n '
            'columns: "chrom", "start", "end", "max_extension" are needed'
        )

    resolution = clr.binsize
    chrom_list = []
    upstream_signal_noise_ratio_list = []
    downstream_signal_noise_ratio_list = []
    ave_bkg_signal_noise_ratio_list = []

    for index, row in regions.iterrows():

        chrom = row['chrom'][3:]
        init_bin = (row['end'] + row['start']) // resolution // 2
        # 20230411 revise
        extension_length = row['max_extension'] * 1000

        if np.isnan(extension_length) or \
                offset >= extension_length:
            tkg_plumb, up_bkg, down_bkg = np.nan, np.nan, np.nan

        else:

            if not chrom in chrom_list:
                chrom_list.append(chrom)
                mat = clr.matrix(balance = norm).fetch(chrom)

            tkg_plumb, up_bkg, down_bkg = set_sampling_box(
                mat=mat, half_width=half_width,
                extension_length=extension_length, offset=offset,
                init_bin=init_bin, resolution=resolution, coverage_ratio=coverage_ratio
            )

        # calculate median interaction in fountain and its signal-over-noise
        # for central area and background areas
        if np.all(
            [isinstance(i, np.ndarray) for i in [up_bkg, down_bkg, tkg_plumb]]
        ):
            bkg_ave = np.nanmean(
                [up_bkg, down_bkg], axis=0
            )

            tkg_plumb_sum = np.nansum(tkg_plumb)
            bkg_ave_sum = np.nansum(bkg_ave)
            bkg_upstream_sum = np.nansum(up_bkg)
            bkg_downstream_sum = np.nansum(down_bkg)

            # calculate signal-noise ratio for average background
            signal_background_ratio = tkg_plumb_sum / bkg_ave_sum

            # calculate signal-noise ratio for upstream background
            signal_upstream_background_ratio = tkg_plumb_sum / bkg_upstream_sum

            # calculate signal-noise ratio for downstream background
            signal_downstream_background_ratio = tkg_plumb_sum / bkg_downstream_sum

            # ratio should not be inf
            if np.isinf(signal_background_ratio):
                signal_background_ratio = np.nan

            if np.isinf(signal_upstream_background_ratio):
                signal_upstream_background_ratio = np.nan

            if np.isinf(signal_downstream_background_ratio):
                signal_downstream_background_ratio = np.nan


            # calculate median of interactions of fountain
            #tkg_plumb_median = np.nanmedian(tkg_plumb_sum)
            # revised 20230411
            tkg_plumb_median = np.nanmedian(tkg_plumb)

        else:
            tkg_plumb_sum = np.nan
            signal_background_ratio = np.nan
            signal_upstream_background_ratio = np.nan
            signal_downstream_background_ratio = np.nan
            tkg_plumb_median = np.nan

        upstream_signal_noise_ratio_list.append(
            signal_upstream_background_ratio
        )

        downstream_signal_noise_ratio_list.append(
            signal_downstream_background_ratio
        )

        ave_bkg_signal_noise_ratio_list.append(
            signal_background_ratio
        )


    regions['signal_noise_upstream'] = upstream_signal_noise_ratio_list
    regions['signal_noise_downstream'] = downstream_signal_noise_ratio_list
    regions['signal_noise_average_background'] = ave_bkg_signal_noise_ratio_list


    return regions

def plumb(
    clr, half_width, extension_length, norm,
    regions, bin_array, offset, coverage_ratio,
    interval_length = 50000, threshold = 0.5
):

    regions = regions.copy()
    regions = regions.sort_values(by = 'chrom').reset_index(drop = True)

    if not set(['chrom', 'start', 'end']).issubset(regions.columns):
        raise TypeError('Invalid dataframe')

    resolution = clr.binsize
    chrom_list = []
    perc_res_list = []
    max_ext_list = []

    for index, row in regions.iterrows():
        chrom = row['chrom'][3:]
        init_bin = (row['end'] + row['start']) // resolution // 2

        if not chrom in chrom_list:
            chrom_list.append(chrom)
            mat = clr.matrix(balance=norm).fetch(chrom)

        max_ext, perc_list = center_background_plumb(
            mat = mat, half_width=half_width,
            extension_length=extension_length, init_bin=init_bin,
            offset = offset, coverage_ratio = coverage_ratio, resolution=resolution,
            bin_array = bin_array, interval_length=interval_length, threshold = threshold
        )

        try:
            perc_list = list(perc_list)
        except:
            perc_list = list(np.zeros_like(bin_array))

        perc_res_list.append(perc_list)

        if not np.isnan(max_ext):

        # append the maximum extension(Kb) to dataframe
            max_ext_length = max_ext // 2 * resolution / 1000

        else:
            max_ext_length = np.nan

        max_ext_list.append(max_ext_length)

    regions['perc_res_list'] = perc_res_list
    regions['max_extension'] = max_ext_list

    return regions


def background_evaluation(
    clr, background_clr, feature, extension_length,
    half_width, offset, norm = 'VC_SQRT'
):
    """
    Compare interactions between Repli-HiC and BL-HiC and
    evaluate the quality of summits. We only leave summits with
    positive value after matrix subtraction.

    Mention: your dataframe has been sorted by chromosome names.
    """
    feature = feature.copy()
    columns = ['chrom', 'start', 'end', 'name', 'SoN', 'strand']

    if not feature.columns.isin(columns).any():
        feature.columns = columns

    feature = feature.sort_values(by = 'chrom').reset_index(drop = True)
    # get resolution
    assert clr.binsize == background_clr.binsize, 'Unmatched resolution'
    resolution = clr.binsize
    half_width = half_width // resolution

    # remove 'chr' label
    if all(['chr' not in i for i in clr.chromnames]):
        if feature.chrom.str.contains('chr').all():
            feature.loc[:, 'chrom'] = [i[3:] for i in feature.chrom]

    chrom_list = []
    val_list = []
    pos_ratio_list = []
    neg_ratio_list = []
    for index, row in feature.iterrows():
        chrom = row['chrom']
        init_bin = (row['end'] + row['start']) // resolution // 2

        if not chrom in chrom_list:
            chrom_list.append(chrom)
            mat = clr.matrix(balance=norm).fetch(chrom)
            background_mat = background_clr.matrix(balance=norm).fetch(chrom)
            mat_subtraction = mat - background_mat

        try:
            extension_length = row['max_extension'] * 1000
        except Exception as e:
            extension_length = extension_length

        target_plumb_val = np.asarray(
            [i for i in plumb_sum(
                mat = mat_subtraction,
                half_width = half_width,
                init_bin = init_bin, extension_length = extension_length,
                resolution = resolution, offset = offset
            )]
        )

        positive_value = target_plumb_val[target_plumb_val>0]
        positive_value_ratio = len(positive_value) / len(target_plumb_val)

        negative_value = target_plumb_val[target_plumb_val<=0]
        negative_value_ratio = len(negative_value) / len(target_plumb_val)


        if len(positive_value) > len(negative_value):
            val_list.append(True)
        else:
            val_list.append(False)

        pos_ratio_list.append(positive_value_ratio)
        neg_ratio_list.append(negative_value_ratio)

    feature.loc[:, 'pos_ratio'] = pos_ratio_list
    feature.loc[:, 'neg_ratio'] = neg_ratio_list
    feature.loc[:, 'subtraction_value'] = val_list

    return feature

