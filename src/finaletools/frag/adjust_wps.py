from __future__ import annotations
from sys import stdin, stdout, stderr
from typing import TextIO, Union
from multiprocessing import Pool
from time import time
import traceback

import numpy as np
import gzip
from numpy.typing import NDArray
import pyBigWig as pbw
from scipy.signal import savgol_filter
from finaletools.utils.utils import _parse_chrom_sizes

def _median_filter(positions: NDArray, data: NDArray, window_size: int):
    running_median = np.array(
        [np.median(data[i:i+window_size]) for i in range(len(data) - window_size)]
    )
    adjusted_data = data[window_size//2 : -(window_size//2)] - running_median
    adjusted_positions = positions[window_size//2 : -(window_size//2)]
    return adjusted_positions, adjusted_data

def _mean_filter(positions: NDArray, data: NDArray, window_size: int):
    running_mean = np.array(
        [np.mean(data[i:i+window_size]) for i in range(len(data) - window_size)]
    )
    adjusted_data = data[window_size//2 : -(window_size//2)] - running_mean
    adjusted_positions = positions[window_size//2 : -(window_size//2)]
    return adjusted_positions, adjusted_data

def _single_adjust_wps_star(args):
    return _single_adjust_wps(*args)


def _single_adjust_wps(
        input_file: str,
        contig: str,
        start: int,
        stop: int,
        median_window_size: int,
        savgol_window_size: int,
        savgol_poly_deg: int,
        mean: bool,
        subtract_edges: bool,
        edge_size: int
):
    try:
        if input_file.endswith('.bw'):
            raw_wps = pbw.open(input_file, 'r')
        else:
            raise ValueError('Invalid filetype for input_file.')
        intervals = np.array(
            list(raw_wps.intervals(contig, start, stop)),
            dtype=[
                ('starts', '<i8'),
                ('stops', '<i8'),
                ('scores', '<f8'),
            ]
        )

        if not all(
            (pos1 + 1 == pos2
             for pos1, pos2
             in zip(intervals['starts'][:-1], intervals['starts'][1:])),
        ):
            raise ValueError(
                'BigWig was found to be nonsequential. There may be '
                'multiple entries for one position or gaps in the '
                'regions specified in the interval file.'
            )

        if subtract_edges:
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
        adjusted_positions = np.zeros((0,), dtype=np.int64)
        stops = np.zeros((0,), dtype=np.int64)
        filtered_scores = np.zeros((0,), dtype=np.float64)

    finally:
        raw_wps.close()

    return (len(adjusted_positions)*[contig],
            adjusted_positions,
            stops,
            filtered_scores)

def adjust_wps(
    input_file: str,
    interval_file: str,
    output_file: str,
    chrom_sizes: str,
    median_window_size: int,
    savgol_window_size: int,
    savgol_poly_deg: int,
    mean: bool,
    subtract_edges: bool,
    edge_size: int,
    workers: int,
    verbose: bool
):
    if verbose:
        start_time = time()
        stderr.write('Reading intervals from bed...\n')

    if interval_file.endswith('.bed') or interval_file.endswith('.bed.gz'):
        end_decrease = median_window_size//2
        intervals = []
        with gzip.open(interval_file, 'rt') if interval_file.endswith('.gz') else open(interval_file, 'rt') as file:
            for line in file:
                contents = line.split('\t')
                contig = contents[0].strip()
                start = int(contents[1])
                stop = int(contents[2])

                if (len(intervals) > 0
                    and intervals[-1][1] == contig
                    and intervals[-1][3] - end_decrease
                    > start + end_decrease
                ):

                    start = intervals[-1][2] 
                    intervals.pop(-1)

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
                    edge_size
                ))
    else:
        raise ValueError('Invalid filetype for the interval file.')

    if verbose:
        stderr.write('Opening pool... \n')
    try:
        pool = Pool(workers)
        processed_scores = pool.imap(_single_adjust_wps_star, intervals)

        if verbose:
            stderr.write('Writing to output\n')

        with pbw.open(output_file, 'w') as output_bw:
            output_bw.addHeader(_parse_chrom_sizes(chrom_sizes))
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
