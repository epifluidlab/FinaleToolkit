from __future__ import annotations
import argparse
import gzip
import time
import os
import tempfile as tf
from multiprocessing.pool import Pool
from typing import Union, TextIO, BinaryIO

import pysam
import numpy as np
from numba import jit
from tqdm import tqdm
from finaletools.utils import (
    frag_bam_to_bed,
    frag_array,
    not_read1_or_low_quality
)

@jit(nopython=True)
def _single_wps(window_start: int,
                window_stop: int,
                window_position: int,
                frag_ends: np.ndarray
                ) -> tuple:
    # count number of totally spanning fragments
    is_spanning = ((frag_ends[:, 0] < window_start)
                   * (frag_ends[:, 1] > window_stop))
    num_spanning = np.sum(is_spanning)

    # count number of fragments with end in window
    is_start_in = ((frag_ends[:, 0] >= window_start)
                   * (frag_ends[:, 0] <= window_stop))
    is_stop_in = ((frag_ends[:, 1] >= window_start)
                  * (frag_ends[:, 1] <= window_stop))
    is_end_in = np.logical_or(is_start_in, is_stop_in)
    num_end_in = np.sum(is_end_in)

    # calculate wps and return
    return (window_position, num_spanning - num_end_in)


@jit(nopython=True)
def _wps_loop(frag_ends: np.ndarray,
              start: int,
              stop: int,
              window_size: int):
    # array to store positions and scores
    scores = np.zeros((stop-start, 2))
    window_centers = np.arange(start, stop, dtype=np.int64)
    scores[:, 0] = window_centers
    window_starts = np.zeros(stop-start)
    window_stops = np.zeros(stop-start)
    np.rint(window_centers - window_size * 0.5, window_starts)
    np.rint(window_centers + window_size * 0.5 - 1, window_stops)
    # inclusive

    for i in range(stop-start):
        scores[i, :] = _single_wps(
            window_starts[i],
            window_stops[i],
            window_centers[i],
            frag_ends)

    return scores


def wps(input_file: Union[str, pysam.AlignmentFile],
        contig: str,
        start: Union[int, str],
        stop: Union[int, str],
        output_file: str=None,
        window_size: int=120,
        fraction_low: int=120,
        fraction_high: int=180,
        quality_threshold: int=30,
        verbose: Union[bool, int]=0
        ) -> np.ndarray:
    """
    Return Windowed Protection Scores as specified in Snyder et al
    (2016) over a region [start,stop).

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM or SAM file containing paired-end fragment reads or its
        path. `AlignmentFile` must be opened in read mode.
    contig : str
    start : int
    stop : int
    output_file : string, optional
    window_size : int, optional
        Size of window to calculate WPS. Default is k = 120, equivalent
        to L-WPS.
    fraction_low : int, optional
        Specifies lowest fragment length included in calculation.
        Default is 120, equivalent to long fraction.
    fraction_high : int, optional
        Specifies highest fragment length included in calculation.
        Default is 120, equivalent to long fraction.
    quality_threshold : int, optional
    workers : int, optional
    verbose : bool, optional

    Returns
    -------
    scores : numpy.ndarray
        np array of shape (n, 2) where column 1 is the coordinate and
        column 2 is the score and n is the number of coordinates in
        region [start,stop)
    """

    if (verbose):
        start_time = time.time()
        print("Reading fragments")

    # set start and stop to ints
    start = int(start)
    stop = int(stop)

    # set minimum and maximum values for fragments. These extend farther
    # than needed
    minimum = round(start - fraction_high)
    maximum = round(stop + fraction_high)

    # read fragments from file
    frag_ends = frag_array(input_file,
                           contig,
                           quality_threshold,
                           minimum=minimum,
                           maximum=maximum,
                           fraction_low=fraction_low,
                           fraction_high=fraction_high,
                           verbose=(verbose>=2))

    if (verbose):
        print("Done reading fragments, preparing for WPS calculation.")
    # check if no fragments exist on this interval
    if (frag_ends.shape == (0, 2)):

        scores = np.zeros((stop-start, 2))
        scores[:, 0] = np.arange(start, stop, dtype=int)
    else:
        scores = _wps_loop(frag_ends, start, stop, window_size)


    # TODO: consider switch-case statements and determine if they
    # shouldn't be used for backwards compatability
    if (type(output_file) == str):   # check if output specified

        if (verbose):
            print('Writing to output file.')

        if output_file.endswith(".wig.gz"): # zipped wiggle
            with gzip.open(output_file, 'wt') as out:
                # declaration line
                out.write(
                    f'fixedStep\tchrom={contig}\tstart={start}\t'
                    f'step={1}\tspan={stop-start}\n'
                    )
                for score in scores[:, 1]:
                    out.write(f'{score}\n')

        elif output_file.endswith(".wig"):  # wiggle
            with open(output_file, 'wt') as out:
                # declaration line
                out.write(
                    f'fixedStep\tchrom={contig}\tstart={start}\tstep='
                    f'{1}\tspan={stop-start}\n'
                    )
                for score in scores[:, 1]:
                    out.write(f'{score}\n')

        else:   # unaccepted file type
            raise ValueError(
                'output_file can only have suffixes .wig or .wig.gz.'
                )

    elif (output_file is not None):
        raise TypeError(
            f'output_file is unsupported type "{type(input_file)}". '
            'output_file should be a string specifying the path of the file '
            'to output scores to.'
            )

    if (verbose):
        end_time = time.time()
        print(f'wps took {end_time - start_time} s to complete')

    return scores