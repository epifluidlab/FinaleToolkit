from __future__ import annotations
from sys import stderr
from typing import Union
import time
import gzip
from os import PathLike

import numpy as np
import pyBigWig as pbw

from finaletoolkit.utils.utils import _get_intervals

def agg_bw(
    input_file: Union[str, PathLike],
    interval_file: Union[str, PathLike],
    output_file: Union[str, PathLike],
    median_window_size: int=1,
    mean: bool=False,
    verbose: bool=False
):
    """
    Takes a BigWig and an interval BED and
    aggregates signal along the intervals with a median filter.

    For aggregating WPS signals, note that the median filter trims the
    ends of each interval by half of the window size of the filter
    while adjusting data. There are two way this can be approached in
    aggregation:

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
        BED file containing intervals. 6th column should have strand.
    output_file : str
    median_window_size : int, optional
        default is 1 (no smoothing). Set to 120 if replicating Snyder et al.
    mean : bool
        use mean filter instead
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
    if (str(interval_file).endswith('.bed')
        or str(interval_file).endswith('.bed.gz')):
        intervals = []
        with (gzip.open(interval_file, 'rt')
              if str(interval_file).endswith('.gz')
              else open(interval_file, 'rt')) as file:
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

    with pbw.open(str(input_file), 'r') as raw_wps: # Path not supported by pbw
        # get size of interval based on first entry in interval_file
        interval_size = intervals[0][2] - intervals[0][1] - median_window_size
        agg_scores = np.zeros(interval_size, dtype=np.int64)
        num_intervals_added = 0
        for contig, start, stop, strand in intervals:
            try:
                signal = raw_wps.values(contig, start, stop)
                if signal is None:
                    print(
                        "There was no information found in the interval: ",
                        contig, start, stop)
                    continue
                values = np.nan_to_num(np.array(signal), nan=0)
            except RuntimeError as e:
                print(e)
                continue

            # trimmed from median filter
            trimmed = values[median_window_size//2:-median_window_size//2]
            if trimmed.shape[0] != interval_size:
                print(
                    f"Trimmed size {trimmed.shape[0]} for {contig}:{start}"
                    f"-{stop} is not equal to "
                    f"interval size {interval_size}. Skipping.")
                continue

            # flip scores if on reverse strand
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

    if str(output_file).endswith('wig'):
        with open(output_file, 'wt') as out:
            if (verbose):
                stderr.write('File opened! Writing...\n')
            # declaration line
            out.write(
                f'fixedStep\tchrom=.\tstart={-interval_size//2}\tstep={1}\t'
                f'span={interval_size}\n'
            )
            for score in agg_scores:
                out.write(f'{score}\n')
    else:
        raise ValueError(
            'The output_file is an unaccepted type. Must be a wiggle file '
            'ending in .wig'
        )
    
    if verbose:
        end_time = time.time()
        stderr.write(f'Aggregating bigWig took {end_time-start_time} s '
                     'to run.\n')

    return agg_scores



