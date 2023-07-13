from __future__ import annotations
from sys import stdin, stdout, stderr
from typing import TextIO, Union
from multiprocessing import Pool
from time import time

import numpy as np
import pyBigWig as pbw

from scipy.signal import savgol_filter

def _median_filter(data, window_size):
    """locally adjusted running median"""
    # Calculate the running median
    running_median = np.array([
        np.median(data[i:i+window_size]) for i in range(len(data) - window_size)
    ])

    # Adjust the data by subtracting the running median
    adjusted_data = data[window_size//2 : -(window_size//2)] - running_median

    return adjusted_data

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
    Golay filter (Savitsky and Golay, 1964) on it.
    """

    """
    adjusted = _median_filter(wps, median_window_size)
    filtered = savgol_filter(adjusted, savgol_window_size, savgol_poly_deg)
    length = len(filtered)
    filt_start = start + median_window_size//2
    """
    # read wps
    try:
        # check input types
        if input_file.endswith('.bw'):
            raw_wps = pbw.open(input_file, 'r')
        else:
            raise ValueError('Invalid filetype for input_file.')

    finally:
        raw_wps.close()

    return None

def process_wps(
    input_file: str,
    interval_file: str,
    output_file: str,
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
    # correct overlaps accounting for median filter
    for interval1, interval2 in intervals[:-1], intervals[1:]:
        if (
            interval1[0] == interval2[0]
            and interval1[2] - end_decrease
            > interval2[1] + end_decrease
        ):
            interval1[2] = interval2[1] + median_window_size

    try:
        pool = Pool(workers)
        processed_scores = pool.imap(_single_process_wps, intervals)
    finally:
        pool.close()








