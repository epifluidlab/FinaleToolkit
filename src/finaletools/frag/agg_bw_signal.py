from __future__ import annotations
from sys import stdin, stdout, stderr
from typing import Union
import time

import numpy as np
import pyBigWig as pbw

def convert_to_array(data):
    data_list=[]
    for element in data:
        data_list.append(element)
    return np.nan_to_num(np.asarray(data_list), nan=0)

def agg_bw_signal(
    input_file: str,
    interval_file: str,
    output_file: str,
    median_window_size: int,
    mean: bool,
    strand_location: int,
    verbose: bool
):
    if verbose:
        start_time = time.time()
        stderr.write('Reading intervals from BED...\n')

    if interval_file.endswith('.bed') or interval_file.endswith('.bed.gz'):
        intervals = []
        with gzip.open(interval_file, 'rt') if interval_file.endswith('.gz') else open(interval_file, 'rt') as file:
            for line in file:
                contents = line.split('\t')
                contig = contents[0]
                start = int(contents[1])
                stop = int(contents[2])
                strand = contents[strand_location]

                intervals.append((
                    contig,
                    int(start),
                    int(stop),
                    strand.strip(),
                ))
    else:
        raise ValueError('Invalid filetype for interval file.')

    with pbw.open(input_file, 'r') as raw_wps:
        interval_size = intervals[0][2] - intervals[0][1] - median_window_size
        agg_scores = np.zeros(interval_size, dtype=np.int64)
        num_intervals_added=0
        for contig, start, stop, strand in intervals:
            try:
                signal = raw_wps.values(contig, start, stop)
                if signal==None:
                    print("There was no information found in the interval: ", contig, start, stop)
                    continue
                values = convert_to_array(signal)
            except RuntimeError as e:
                print(e)
                continue
            if median_window_size == 0:
                trimmed = values
            else:
                trimmed = values[median_window_size//2:-median_window_size//2]
            if trimmed.shape[0] != interval_size:
                print(f"Trimmed size {trimmed.shape[0]} is not equal to interval size {interval_size}. Skipping.")
                continue
            if strand == '+':
                agg_scores = agg_scores + trimmed
                num_intervals_added+=1
            elif strand == '-':
                agg_scores = agg_scores + np.flip(trimmed)
                num_intervals_added+=1
            elif verbose:
                stderr.write(
                    'A segment without strand was encountered. Skipping.'
                )
    if mean:
        agg_scores = agg_scores/num_intervals_added
    positions = np.arange(-interval_size//2, interval_size//2)
    if output_file.endswith('wig'):
        with open(output_file, 'wt') as out:
            if (verbose):
                stderr.write(f'File opened! Writing...\n')
            out.write(
                f'fixedStep\tchrom=.\tstart={-interval_size//2}\tstep={1}\tspan'
                f'={interval_size}\n'
            )
            for score in agg_scores:
                out.write(f'{score}\n')
    else:
        raise ValueError(
            'The output_file is an unaccepted type. Must be a wiggle file ending in .wig'
        )
    
    if verbose:
        end_time = time.time()
        stderr.write(f'Aggregating bigWig signals took {end_time-start_time} s to run.\n')

    return agg_scores



