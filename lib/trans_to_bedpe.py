import pandas as pd

def dataframe_to_bedpe(regions, resolution):

    up_coords_s, down_coords_s = [], []
    up_coords_e, down_coords_e = [], []
    chrom_list = []
    for index, row in regions.iterrows():
        chrom = row['chrom']

        mid_coord_s = int(row['start'] + row['end']) // 2

        up_coord_s = int(mid_coord_s - row['max_extension'] * 1000)
        up_coord_e = int(up_coord_s + resolution)

        down_coord_s = int(mid_coord_s + row['max_extension'] * 1000)
        down_coord_e = int(down_coord_s + resolution)

        chrom_list.append(chrom)

        up_coords_s.append(up_coord_s)
        up_coords_e.append(up_coord_e)

        down_coords_s.append(down_coord_s)
        down_coords_e.append(down_coord_e)

    df = pd.DataFrame(
        {'chr1': chrom_list,
         'x1': up_coords_s,
         'x2': up_coords_e,
         'chr2': chrom_list,
         'y1': down_coords_s,
         'y2': down_coords_e}
    )

    return df


def label_sampling_box(regions, resolution, half_width):

    center_upstream_bound_list, center_downstream_bound_list = [], []
    upstream_upper_bound_list, downstream_lower_bound_list = [], []
    chrom_list = []
    for index, row in regions.iterrows():

        center = row['start']
        chrom = row['chrom']

        center_upstream_bound = center - half_width * resolution
        center_downstream_bound = center + half_width * resolution

        upstream_upper_bound = center_upstream_bound - 2 * half_width * resolution
        downstream_lower_bound = center_downstream_bound + 2 * half_width * resolution

        center_upstream_bound_list.append(center_upstream_bound)
        center_downstream_bound_list.append(center_downstream_bound)
        upstream_upper_bound_list.append(upstream_upper_bound)
        downstream_lower_bound_list.append(downstream_lower_bound)
        chrom_list.append(chrom)

    center_upstream_bound_list = np.asarray(center_upstream_bound_list)
    center_downstream_bound_list = np.asarray(center_downstream_bound_list)
    upstream_upper_bound_list = np.asarray(upstream_upper_bound_list)
    downstream_lower_bound_list = np.asarray(downstream_lower_bound_list)

    df_upstream_upper_bound_list = pd.DataFrame(
        {'chr1': chrom_list,
         'x1': upstream_upper_bound_list,
         'x2': upstream_upper_bound_list + resolution,
         'chr2': chrom_list,
         'y1': upstream_upper_bound_list,
         'y2': upstream_upper_bound_list + resolution}
    )

    df_downstream_lower_bound_list = pd.DataFrame(
        {'chr1': chrom_list,
         'x1': downstream_lower_bound_list,
         'x2': downstream_lower_bound_list + resolution,
         'chr2': chrom_list,
         'y1': downstream_lower_bound_list,
         'y2': downstream_lower_bound_list + resolution}
    )

    df_center_upstream_bound_list = pd.DataFrame(
        {'chr1': chrom_list,
         'x1': center_upstream_bound_list,
         'x2': center_upstream_bound_list + resolution,
         'chr2': chrom_list,
         'y1': center_upstream_bound_list,
         'y2': center_upstream_bound_list + resolution}
    )

    df_center_downstream_bound_list = pd.DataFrame(
        {'chr1': chrom_list,
         'x1': center_downstream_bound_list,
         'x2': center_downstream_bound_list + resolution,
         'chr2': chrom_list,
         'y1': center_downstream_bound_list,
         'y2': center_downstream_bound_list + resolution}
    )

    df_merged = pd.concat(
        [df_upstream_upper_bound_list, df_downstream_lower_bound_list,
         df_center_upstream_bound_list, df_center_downstream_bound_list]
    )

    return df_merged
