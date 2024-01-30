def generate_summits_bed(clr, track, out_dir):

    resolution = clr.binsize
    suffix = '_' + str(resolution // 1000) + 'kb.bed'

    out_dir = out_dir + 'SoN_summits/'
    os.system('mkdir %s' % out_dir)

    for chrom in clr.chromnames:
        chrom = 'chr' + chrom

        # get SoN signal tracks for targeted chromosome
        signal_filter = track[track[0] == chrom]
        signal_filter_value = signal_filter[3].values

        # When find summits, we ignore the negative SoN values
        signal_filter_value[signal_filter_value < 0] = 0
        signal_filter[3] = signal_filter_value

        # get positions of summits
        poss, proms = find_peak_prominence(signal_filter[3].values)

        chr_cord_start = [i * resolution for i in poss]
        chr_cord_end = [(i + 1) * resolution for i in poss]
        names = np.repeat('.', len(poss))
        strands = np.repeat('.', len(poss))

        # get SoN score
        SoN_values = signal_filter[3].values[poss]

        data_dict = {
            'chr': np.repeat(chrom, len(poss)),
            'start': chr_cord_start,
            'end': chr_cord_end,
            'name': names,
            'score': SoN_values,
            'strand': strands
        }

        output = chrom + suffix
        df = pd.DataFrame(data_dict)
        df.to_csv(
            out_dir + output,
            sep='\t',
            header=None,
            index=None
        )

    # merge all summits into one file
    os.system(
        'cat %s*.bed > %sSummits_%s_merged.bed' % (out_dir, out_dir, resolution)
    )

    # remove summits for individual chromosome
    os.system(
        'rm %schr*' % out_dir
    )
