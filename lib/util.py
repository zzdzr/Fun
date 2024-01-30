import subprocess
import bioframe
import pandas as pd

def align_track_with_chromsize(track, chromsizes):
    """
    Aligns a signal track with chromosome sizes, ensuring that 'end' coordinates do not exceed chromosome sizes.

    Args:
        track (DataFrame): Signal track of SoN.
        chromsizes (DataFrame): DataFrame containing information of chromosome sizes.

    Returns:
        DataFrame: Refined SoN with adjusted 'end' coordinates.
    """
    track = track.copy()

    # Get chromosome information from the track
    chrom = track['chr'].unique().item()

    # Get chromosome size
    size = chromsizes.at[chrom]

    # Update 'end' coordinates exceeding chromosome size
    track.loc[track['end'] > size, 'end'] = size

    return track


def bedgraph_to_bigwig(bedGraphToBigWig, track_path, chromsizes_path, bigwig_output_path):
    """

    Args:
        bedGraphToBigWig (str): Path to bedGraphToBigWig from UCSC
        bedgraph_path (str): Path to the input bedGraph file.
        chromsize_path (str): Path to the chromosome size file.
        bigwig_output_path (str): Path to save the output BigWig file.

    Returns:
        None
    """
    cmd = [bedGraphToBigWig, track_path, chromsizes_path, bigwig_output_path]
    subprocess.run(cmd)