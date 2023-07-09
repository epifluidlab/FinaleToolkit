from __future__ import annotations
from sys import stdin, stdout, stderr
from typing import TextIO

import numpy as np

from scipy.signal import savgol_filter

def _larm(data, window_size):
    """locally adjusted running median"""
    # Calculate the running median
    running_median = np.array([
        np.median(data[i:i+window_size]) for i in range(len(data) - window_size)
    ])

    # Adjust the data by subtracting the running median
    adjusted_data = data[window_size//2 : -(window_size//2)] - running_median

    return adjusted_data

def process_wps(
        input_file: str,
        output_file: str,
        larm_window_size: int=1000,
        savgol_window_size: int=21,
        savgol_poly_deg: int=2,
    ):
    """
    Takes a wps WIG file and applies a median filter and a Savitsky-
    Golay filter (Savitsky and Golay, 1964) on it.
    """
    try:
        # read from file or stdin
        if input_file == '-':
            wig_file = stdin
        else:
            wig_file = open(input_file, 'r')

        wps = []
        is_first_line = True
        for line in wig_file:
            if line.isspace():
                break
            # Read info from header
            if is_first_line:
                contents = line.split('\t')
                for item in contents:
                    if 'chrom' in item:
                        chrom = item.split('=')[1].strip()
                    if 'start' in item:
                        start = int(item.split('=')[1])
                    if 'span' in item:
                        span = int(item.split('=')[1])

                is_first_line = False
                continue
            else:
                wps.append(float(line))
    finally:
        if input_file != '-':
            wig_file.close()

    adjusted = _larm(wps, larm_window_size)
    filtered = savgol_filter(adjusted, savgol_window_size, savgol_poly_deg)
    length = len(filtered)
    filt_start = start + larm_window_size//2

    # write to file or stdout
    try:
        if output_file == '-':
            wig_file = stdout
        else:
            wig_file = open(output_file, 'w')

        header = (
            f'fixedStep\tchrom={chrom}\tstart={filt_start}\tstep=1\tspan='
            f'{length}\n'
        )
        wig_file.write(header)
        for score in filtered:
            wig_file.write(f'{score}\n')
    finally:
        if output_file != '-':
            wig_file.close()

    


