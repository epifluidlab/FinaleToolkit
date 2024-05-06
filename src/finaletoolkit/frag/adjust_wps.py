from __future__ import annotations
from sys import stdin, stdout, stderr
from typing import TextIO, Union
from multiprocessing import Pool
from time import time
import traceback
import gzip

import numpy as np
from numpy.typing import NDArray
import pyBigWig as pbw
from scipy.signal import savgol_filter

from finaletoolkit.utils import genome2list

def _median_filter(positions: NDArray, data: NDArray, window_size: int):
    """locally adjusted running median"""
    # Calculate the running median
    running_median = np.array(
        [np.median(data[i:i+window_size]) for i in range(len(data) - window_size)]
    )

    # Adjust the data by subtracting the running median
    adjusted_data = data[window_size//2 : -(window_size//2)] - running_median

    # trim positions
    adjusted_positions = positions[window_size//2 : -(window_size//2)]

    return adjusted_positions, adjusted_data


def _mean_filter(positions: NDArray, data: NDArray, window_size: int):
    """locally adjusted running mean"""
    # Calculate the running mean
    running_mean = np.array(
        [np.mean(data[i:i+window_size]) for i in range(len(data) - window_size)]
    )

    # Adjust the data by subtracting the running median
    adjusted_data = data[window_size//2 : -(window_size//2)] - running_mean

    # trim positions
    adjusted_positions = positions[window_size//2 : -(window_size//2)]

    return adjusted_positions, adjusted_data


def _single_adjust_wps_star(args):
    """
    Version of single_process_wps that accepts a tuple of args. Used for
    multiprocess.imap.
    """
    return _single_adjust_wps(*args)


def _single_adjust_wps(
        input_file: str,
        contig: str,
        start: int,
        stop: int,
        median_window_size: int=1000,
        savgol_window_size: int=21,
        savgol_poly_deg: int=2,
        mean: bool=False,
        subtract_edges: bool=False,
        edge_size: int=500,
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

        # scores are read from BigWig into an numpy structured array
        intervals = np.array(
            list(raw_wps.intervals(contig, start, stop)),
            dtype=[
                ('starts', '<i8'),
                ('stops', '<i8'),
                ('scores', '<f8'),
            ]
        )

        # check BigWig for errors
        if not all(
            (pos1 + 1 == pos2
             for pos1, pos2
             in zip(intervals['starts'][:-1], intervals['starts'][1:])),
        ):
            # TODO: create special error for invalid file formats
            raise ValueError(
                'BigWig was found to be nonsequential. There may be '
                'multiple entries for one position or gaps in the '
                'regions specified in the interval file.'
            )

        # adjusting/filtering
        if subtract_edges:
            # TODO: add option for edge size
            start_mean = np.mean(intervals['scores'][:edge_size])
            stop_mean = np.mean(intervals['scores'][-edge_size:])
            mean = np.mean([start_mean, stop_mean])
            intervals['scores'] = intervals['scores'] - mean

        if not mean:
            adjusted_positions, adjusted_scores = _median_filter(
                intervals['starts'], intervals['scores'], median_window_size)
        else:
            adjusted_positions, adjusted_scores = _mean_filter(
                intervals['starts'], intervals['scores'], median_window_size)

        filtered_scores = savgol_filter(
            adjusted_scores, savgol_window_size, savgol_poly_deg)

        assert len(adjusted_positions) == len(filtered_scores)

        stops = adjusted_positions+1

    except RuntimeError as e:
        traceback.print_exception(e)
        stderr.write(
            f'Invalid interval detected:\n'
            f'{contig}:{start}-{stop}. This interval will be skipped.\n'
        )
        # While writing to output, empty intervals like this are skipped
        adjusted_positions = np.zeros((0,), dtype=np.int64)
        stops = np.zeros((0,), dtype=np.int64)
        filtered_scores = np.zeros((0,), dtype=np.float64)

    finally:
        raw_wps.close()

    # the contig is stored in a list while the rest of the values are
    # in NDArrays. This may change in future versions.
    return (len(adjusted_positions)*[contig],
            adjusted_positions,
            stops,
            filtered_scores)

def adjust_wps(
    input_file: str,
    interval_file: str,
    output_file: str,
    genome_file: str,
    median_window_size: int=1000,
    savgol_window_size: int=21,
    savgol_poly_deg: int=2,
    mean: bool=False,
    subtract_edges: bool=False,
    edge_size: int=500,
    workers: int=1,
    verbose: Union(bool, int)=False
):
    """
    Adjusts raw WPS data in a BigWig by applying a median filter and
    Savitsky-Golay filter (Savitsky and Golay, 1964).

    Parameters
    ----------
    input_file : str
        Path string to a BigWig containing raw WPS data.
    interval_file : str
        BED format file containing intervals over which WPS was
        calculated on.
    output_file : str
        BigWig file to write adjusted WPS to.
    genome_file : str
        The genome file for the reference genome that WGS was aligned
        to. A tab delimited file where column 1 contains the name of
        chromosomes and column 2 contains chromosome length.
    median_window_size : int, optional
        Size of median filter window. Default is 1000.
    savgol_window_size : int, optional
        Size of Savitsky Golay filter window. Default is 21.
    savgol_poly_deg : int, optional
        Degree polynomial for Savitsky Golay filter. Default is 2.
    mean : bool, optional
        If true, a mean filter is used instead of median. Default is
        False.
    subtract_edges : bool, optional
        If true, take the median of the first and last 500 bases in a
        window and subtract from the whole interval. Default is False.
    workers : int, optional
        Number of processes to use. Default is 1.
    verbose : bool or int, optional
        Default is False.

    """
    if verbose:
        start_time = time()
        stderr.write('Reading intervals from bed...\n')

    # read intervals
    if interval_file.endswith('.bed') or interval_file.endswith('.bed.gz'):
        # amount taken by median filter
        end_decrease = median_window_size//2
        intervals = []
        with (gzip.open(interval_file, 'rt')
              if interval_file.endswith('.gz')
              else open(interval_file, 'rt')) as file:
            for line in file:
                # read segment from BED
                contents = line.split('\t')
                contig = contents[0].strip()
                start = int(contents[1])
                stop = int(contents[2])

                # Checks for overlap with previous interval, accounting
                # for change in interval size from median filter.
                # This is needed to avoid duplicate entries in the
                # BigWig file.
                if (len(intervals) > 0
                    and intervals[-1][1] == contig
                    and intervals[-1][3] - end_decrease
                    > start + end_decrease
                ):
                    # if overlap, pop first from list and append the
                    # union of the intervals to list
                    start = intervals[-1][2]    # start of first
                    intervals.pop(-1)   # pop first

                intervals.append((
                    input_file,
                    contig,
                    int(start),
                    int(stop),
                    median_window_size,
                    savgol_window_size,
                    savgol_poly_deg,
                    mean,
                    subtract_edges,
                    edge_size,
                ))
    else:
        raise ValueError('Invalid filetype for interval_file.')

    if verbose:
        stderr.write('Opening pool...\n')
        
    try:
        # use pool of processes to process wps scores into an iterator
        pool = Pool(workers)
        processed_scores = pool.imap(_single_adjust_wps_star, intervals)

        if verbose:
            stderr.write('Writing to output\n')

        # write to output
        with pbw.open(output_file, 'w') as output_bw:
            output_bw.addHeader(genome2list(genome_file))
            for scores in processed_scores:
                contigs, starts, stops, values = scores
                if len(contigs) == 0:
                    continue
                try:
                    output_bw.addEntries(
                        contigs,
                        starts,
                        ends=stops,
                        values=values,
                    )
                except RuntimeError as e:
                    traceback.print_exception(e)
                    stderr.write(
                        f'RuntimeError encountered while writing to '
                        f'{output_file} at interval {contigs[[0]]}:'
                        f'{starts[0]}-{stops[-1]}\n'
                    )

    finally:
        pool.close()

    if verbose:
        end_time = time()
        stderr.write(f'Adjust-WPS took {end_time-start_time} s to run.\n')








