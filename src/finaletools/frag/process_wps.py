from __future__ import annotations
from sys import stdin, stdout, stderr
from typing import TextIO, Union
from multiprocessing import Pool
from time import time

import numpy as np
from numpy.typing import NDArray
import pyBigWig as pbw
from scipy.signal import savgol_filter

from finaletools.utils import genome2list

def _median_filter(positions: NDArray, data: NDArray, window_size: int):
    """locally adjusted running median"""
    # Calculate the running median
    running_median = np.array([
        np.median(data[i:i+window_size]) for i in range(len(data) - window_size)
    ])

    # Adjust the data by subtracting the running median
    adjusted_data = data[window_size//2 : -(window_size//2)] - running_median

    # trim positions
    adjusted_positions = positions[window_size//2 : -(window_size//2)]

    return adjusted_positions, adjusted_data


def _single_process_wps_star(args):
    """
    Version of single_process_wps that accepts a tuple of args. Used for
    multiprocess.imap.
    """
    return _single_process_wps(*args)


def _single_process_wps(
        input_file: str,
        contig: str,
        start: int,
        stop: int,
        median_window_size: int=1000,
        savgol_window_size: int=21,
        savgol_poly_deg: int=2,
):
    """
    Takes a wps WIG file and applies a median filter and a Savitsky-
    Golay filter (Savitsky and Golay, 1964) on it. Also returns
    positions for the window (which is truncated by the median filter).
    """

    # read scores
    try:
        # check input types
        if input_file.endswith('.bw'):
            raw_wps = pbw.open(input_file, 'r')
        else:
            raise ValueError('Invalid filetype for input_file.')

        intervals = np.array(
            raw_wps.intervals(contig, start, stop),
            dtype=[
                ('starts', '<i8'),
                ('stops', '<i8'),
                ('scores', '<f8'),
            ]
        )

        if not all(
            (pos1 + 1 == pos2
             for pos1, pos2
             in zip(intervals['starts',:-1], intervals['starts',1:])),
        ):
            # TODO: create special error for invalid file formats
            raise ValueError(
                'BigWig was found to be nonsequential. There may be '
                'multiple entries for one position or gaps in the '
                'regions specified in the interval file.'
            )

        adjusted_positions, adjusted_scores = _median_filter(
            intervals['starts'], intervals['scores'], median_window_size)

        filtered_scores = savgol_filter(
            adjusted_scores, savgol_window_size, savgol_poly_deg)

    finally:
        raw_wps.close()

    assert len(adjusted_positions) == len(filtered_scores)

    stops = adjusted_positions+1

    # the contig is stored in a list while the rest of the values are
    # in NDArrays. This may change in future versions.
    return (len(adjusted_positions)*[contig],
            adjusted_positions,
            stops,
            filtered_scores)

def process_wps(
    input_file: str,
    interval_file: str,
    output_file: str,
    genome_file: str,
    median_window_size: int=1000,
    savgol_window_size: int=21,
    savgol_poly_deg: int=2,
    workers: int=1,
    verbose: Union(bool, int)=False
):
    if verbose:
        start_time = time()

    # read intervals
        # read intervals
    if interval_file.endswith('.bed') or interval_file.endswith('.bed.gz'):
        intervals = []
        with open(interval_file, 'r') as file:
            for line in file:
                contents = line.split('\t')
                contig = contents[0]
                start = int(contents[1])
                stop = int(contents[2])
                intervals.append(
                    input_file,
                    contig,
                    start,
                    stop,
                    median_window_size,
                    savgol_window_size,
                    savgol_poly_deg
                )
    else:
        raise ValueError('Invalid filetype for interval_file.')

    # amount taken by median filter
    end_decrease = median_window_size//2
    # correct overlaps accounting for shortening from median filter
    for interval1, interval2 in intervals[:-1], intervals[1:]:
        if (
            interval1[0] == interval2[0]    # same contig
            and interval1[2] - end_decrease
            > interval2[1] + end_decrease
        ):
            interval1[2] = interval2[1] + median_window_size

    try:
        # use pool of processes to process wps scores into an iterator
        pool = Pool(workers)
        processed_scores = pool.imap(_single_process_wps_star, intervals)

        # read to output
        with pbw.open(output_file, 'w') as output_bw:
            output_bw.addHeader(genome2list(genome_file))
            for scores in processed_scores:
                contigs, starts, stops, values = processed_scores
                if len(contigs) == 0:
                    continue
                try:
                    output_bw.addEntries(
                        contigs,
                        starts,
                        ends=stops,
                        values=values,
                    )
                except RuntimeError:
                    stderr.write(
                        f'RuntimeError encountered while writing to '
                        f'{output_file} at interval {contigs[[0]]}:'
                        f'{starts[0]}-{stops[-1]}\n'
                    )

    finally:
        pool.close()

    if verbose:
        end_time = time()
        stderr.write(f'Process-WPS took {end_time-start_time} s to run.\n')








