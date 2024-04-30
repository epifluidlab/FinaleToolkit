from __future__ import annotations
from sys import stdin, stdout, stderr
from typing import Union
import time

import numpy as np
import pyBigWig as pbw

def agg_wps(
    input_file: str,
    interval_file: str,
    output_file: str,
    median_window_size: int=0,
    verbose: Union[bool, int]=False
):
    """
    Takes a BigWig containing adjusted WPS scores and an interval BED and
    aggregates scores along the intervals.

    Note that the median filter trims the ends of each interval by half
    of the window size of the filter while adjusting raw WPS data. There
    are two way this can be approached in aggregation:
    1. supply an interval file containing smaller intervals. e.g. if
    you used 5kb intervals for WPS and used a median filter window
    of 1kb, supply a BED file with 4kb windows to this function.

    2. provide the size of the median filter window in
    `median_window_size` along with the original intervals. e.g if
    5kb intervals were used for WPS and a 1kb median filter window
    was used, supply the 5kb bed file and median filter window size
    to this function.

    Do not do both of these at once.

    Parameters
    ----------
    input_file : str
    interval_file : str
    output_file : str
    median_window_size : int, optional
        default is 0
    verbose : int or bool, optional
        default is False

    Return
    ------
    agg_scores : NDArray
    """
    if verbose:
        start_time = time.time()
        stderr.write('Reading intervals from bed...\n')

    # reading intervals from interval_file into a list
    if interval_file.endswith('.bed') or interval_file.endswith('.bed.gz'):
        # TODO: add support for bed.gz
        if interval_file.endswith('.gz'):
            raise NotImplementedError('bed.gz not supported yet')
        intervals = []
        with open(interval_file, 'r') as file:
            for line in file:
                # read segment from BED
                contents = line.split('\t')
                contig = contents[0]
                start = int(contents[1])
                stop = int(contents[2])
                strand = contents[5]

                intervals.append((
                    contig,
                    int(start),
                    int(stop),
                    strand.strip(),
                ))
    else:
        raise ValueError('Invalid filetype for interval_file.')

    with pbw.open(input_file, 'r') as raw_wps:
        # get size of interval based on first entry in interval_file
        interval_size = intervals[0][2] - intervals[0][1] - median_window_size
        agg_scores = np.zeros(interval_size, dtype=np.int64)
        # HACK: excepting for random intervals with wrong size or improperly
        # formatted BigWigs
        for contig, start, stop, strand in intervals:
            try:
                values = np.array(raw_wps.values(contig, start, stop))
            except RuntimeError as e:
                print(e)
                continue
            # trimmed from median filter
            trimmed = values[median_window_size//2:-median_window_size//2]
            if trimmed.shape[0] != interval_size:
                print(f"Trimmed size {trimmed.shape[0]} is not equal to interval size {interval_size}. Skipping.")
                continue
            # flip scores if on reverse strand
            if strand == '+':
                agg_scores = agg_scores + trimmed
            elif strand == '-':
                agg_scores = agg_scores + np.flip(trimmed)
            elif verbose:
                stderr.write(
                    'A segment without strand was encountered. Skipping.'
                )
    positions = np.arange(-interval_size//2, interval_size//2)
    if output_file.endswith('wig'):
        with open(output_file, 'wt') as out:
            if (verbose):
                stderr.write(f'File opened! Writing...\n')
            # declaration line
            out.write(
                f'fixedStep\tchrom=.\tstart={-interval_size//2}\tstep={1}\tspan'
                f'={interval_size}\n'
            )
            for score in agg_scores:
                out.write(f'{score}\n')
    else:
        raise ValueError(
            'output_file is unaccepted type.'
        )
    
    if verbose:
        end_time = time.time()
        stderr.write(f'Agg-WPS took {end_time-start_time} s to run.\n')

    return agg_scores



