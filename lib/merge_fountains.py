import numpy as np

def merge_fountains(regions, max_distance = 50000):

    regions = regions.copy()

    # chrom list
    chroms = regions['chrom'].unique()
    # list of dataframe index for storing the best fountain candidates
    unique_index_list = []

    for chrom in chroms:
        df = regions[regions['chrom'] == chrom]
        df_starts = np.asarray(df['start'], dtype = int)

        upstream_bound = df_starts - max_distance
        downstream_bound = df_starts + max_distance

        for up, down in zip(upstream_bound, downstream_bound):
            df_overlap = df[
                (df['start'] >= up) & (df['start'] <= down)
            ]

            # record the dataframe index for the best fountain candidate
            if len(df_overlap) > 1:
                df_overlap_score = df_overlap['signal_noise_average_background'].values
                df_overlap_row_index = list(df_overlap.index)

                max_idx = np.argmax(df_overlap_score)
                df_max_idx = df_overlap.iloc[max_idx].name

                # if the max index has been stored in the list
                if not df_max_idx in unique_index_list:

                    try:
                        if unique_index_list[-1] in df_overlap_row_index:
                            unique_index_list.pop()
                    except IndexError:
                        unique_index_list.append(df_overlap.index[0])

                    unique_index_list.append(df_max_idx)

            else:
                # if this fountain is unique in this region,
                # just record the index of this line
                unique_index_list.append(df_overlap.index[0])

    return regions.iloc[unique_index_list]