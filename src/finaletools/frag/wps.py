from __future__ import annotations
import argparse
import gzip
import time
import os
import tempfile as tf
from multiprocessing.pool import Pool
from typing import Union, TextIO, BinaryIO
from sys import stdout, stderr, getsizeof

import pysam
import numpy as np
from numba import jit
from tqdm import tqdm

from finaletools.utils.utils import frag_array

@jit(nopython=True)
def _single_wps(contig: str,
                window_start: int,
                window_stop: int,
                window_position: int,
                frag_ends: np.ndarray
                ) -> tuple:
    is_spanning = ((frag_ends["start"] < window_start)
                   * (frag_ends["stop"] > window_stop))
    num_spanning = np.sum(is_spanning)

    is_start_in = ((frag_ends["start"] >= window_start)
                   * (frag_ends["start"] <= window_stop))
    is_stop_in = ((frag_ends["stop"] >= window_start)
                  * (frag_ends["stop"] <= window_stop))
    is_end_in = np.logical_or(is_start_in, is_stop_in)
    num_end_in = np.sum(is_end_in)

    return (contig, window_position, num_spanning - num_end_in)


def wps(input_file: Union[str, pysam.AlignmentFile],
        contig: str,
        start: Union[int, str],
        stop: Union[int, str],
        output_file: str,
        window_size: int,
        min_length: int,
        max_length: int,
        intersect_policy: str,
        quality_threshold: int,
        verbose: bool
        ) -> np.ndarray:

    if (verbose):
        start_time = time.time()
        stderr.write("[finaletools-wps] Reading fragments\n")
        stderr.write(f'Region: {contig}:{start}-{stop}\n')

    start = int(start)
    stop = int(stop)

    minimum = max(round(start - min_length), 0)
    maximum = round(stop + min_length)

    frag_ends = frag_array(input_file,
                           contig,
                           quality_threshold,
                           start=minimum,
                           stop=maximum,
                           min=min_length,
                           max=max_length,
                           intersect_policy=intersect_policy,
                           verbose=(verbose>=2))

    if (verbose):
        stderr.write("Done reading fragments, preparing for WPS calculation.\n")
    
    if (frag_ends.shape == (0)):
        scores = np.zeros(
            stop-start,
            dtype=[
                ('contig', 'U16'),
                ('start', 'i8'),
                ('wps', 'i8'),
            ]
        )
        scores['start'] = np.arange(start, stop, dtype=int)
        scores['contig'] = contig
    else:
        window_centers = np.arange(start, stop, dtype=np.int64)
        window_starts = np.zeros(stop-start)
        window_stops = np.zeros(stop-start)
        np.rint(window_centers - window_size * 0.5, window_starts)
        np.rint(window_centers + window_size * 0.5 - 1, window_stops)

        scores = np.zeros(
            stop-start,
            dtype=[
                ('contig', 'U16'),
                ('start', 'i8'),
                ('wps', 'i8'),
            ]
        )
        for i in range(stop-start):
            scores[i] = _single_wps(
                contig,
                window_starts[i],
                window_stops[i],
                window_centers[i],
                frag_ends
            )

    if (type(output_file) == str): 

        if (verbose):
            stderr.write('Writing to output file.\n')

        if output_file.endswith(".wig.gz"):
            with gzip.open(output_file, 'wt') as out:
                out.write(
                    f'fixedStep\tchrom={contig}\tstart={start}\t'
                    f'step={1}\tspan={stop-start}\n'
                    )
                for score in scores['wps']:
                    out.write(f'{score}\n')

        elif output_file.endswith(".wig"):  # wiggle
            with open(output_file, 'wt') as out:
                out.write(
                    f'fixedStep\tchrom={contig}\tstart={start}\tstep='
                    f'{1}\tspan={stop-start}\n'
                    )
                for score in scores['wps']:
                    out.write(f'{score}\n')

        elif output_file == '-':
            stdout.write(
                f'fixedStep\tchrom={contig}\tstart={start}\tstep='
                f'{1}\tspan={stop-start}\n'
                )
            for score in scores['wps']:
                stdout.write(f'{score}\n')
            stdout.flush()

        else:   # unaccepted file type
            raise ValueError(
                'The output file can only have suffixes .wig or .wig.gz.'
                )

    elif (output_file is not None):
        raise TypeError(
            f'The output file is unsupported type "{type(input_file)}". '
            'Output file should be a string specifying the path of the file '
            'to output scores to.'
            )

    if (verbose):
        end_time = time.time()
        stderr.write(f'WPS took {end_time - start_time} s to complete\n')

    return scores
