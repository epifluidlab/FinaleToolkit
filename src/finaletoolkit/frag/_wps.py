from __future__ import annotations
import gzip
import time
from typing import Union
from sys import stdout, stderr
import warnings

import pysam
import numpy as np
from numba import jit

from finaletoolkit.utils import frag_array
from ..utils.typing import ChromSizes

@jit(nopython=True)
def _single_nt_wps(chrom: str,
                window_start: int,
                window_stop: int,
                window_position: int,
                frag_ends: np.ndarray
                ) -> tuple:
    # count number of totally spanning fragments
    is_spanning = ((frag_ends["start"] < window_start)
                   * (frag_ends["stop"] > window_stop))
    num_spanning = np.sum(is_spanning)

    # count number of fragments with end in window
    is_start_in = ((frag_ends["start"] >= window_start)
                   * (frag_ends["start"] <= window_stop))
    is_stop_in = ((frag_ends["stop"] >= window_start)
                  * (frag_ends["stop"] <= window_stop))
    is_end_in = np.logical_or(is_start_in, is_stop_in)
    num_end_in = np.sum(is_end_in)

    # calculate wps and return
    return (chrom, window_position, num_spanning - num_end_in)


def wps(input_file: Union[str, pysam.AlignmentFile],
        chrom: str,
        start: int,
        stop: int,
        chrom_size: int,
        output_file: str | None = None,
        window_size: int = 120,
        min_length: int = 120,
        max_length: int = 180,
        quality_threshold: int = 30,
        verbose: bool | int = 0,
        fraction_low: int | None = None,
        fraction_high: int | None = None,
        ) -> np.ndarray:
    """
    Return (raw) Windowed Protection Scores as specified in Snyder et al
    (2016) over a region [start,stop).

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, CRAM or tabix file containing paired-end fragment reads or its
        path. `AlignmentFile` must be opened in read mode.
    chrom : str
    start : int
    stop : int
    chrom_size : int
        Size of chrom
    output_file : string, optional
    window_size : int, optional
        Size of window to calculate WPS. Default is k = 120, equivalent
        to L-WPS.
    min_length : int, optional
        Specifies lowest fragment length included in calculation.
        Default is 120, equivalent to long WPS.
    max_length : int, optional
        Specifies highest fragment length included in calculation.
        Default is 180, equivalent to long WPS.
    quality_threshold : int, optional
    workers : int, optional
    verbose : bool, optional
    fraction_low : int, optional
        Deprecated alias for `min_length`
    fraction_high : int, optional
        Deprecated alias for `max_length`

    Returns
    -------
    scores : numpy.ndarray
        np struct array of with columns `contig`, `start`, and `wps`.
    """

    if (verbose):
        start_time = time.time()
        stderr.write("[finaletoolkit-wps] Reading fragments\n")
        stderr.write(f'Region: {chrom}:{start}-{stop}\n')
    
    # Pass aliases and check for conflicts
    if fraction_low is not None and min_length is None:
        min_length = fraction_low
        warnings.warn("fraction_low is deprecated. Use min_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
    elif fraction_low is not None and min_length is not None:
        warnings.warn("fraction_low is deprecated. Use min_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
        raise ValueError(
            'fraction_low and min_length cannot both be specified')

    if fraction_high is not None and max_length is None:
        max_length = fraction_high
        warnings.warn("fraction_high is deprecated. Use max_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
    elif fraction_high is not None and max_length is not None:
        warnings.warn("fraction_high is deprecated. Use max_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
        raise ValueError(
            'fraction_high and max_length cannot both be specified')

    # set start and stop to ints
    start = int(start)
    stop = int(stop)

    # set minimum and maximum values for fragments. These extend farther
    # than needed
    minimum = max(round(start - max_length), 0)
    maximum = min(round(stop + max_length), chrom_size)

    # read fragments from file
    frag_ends = frag_array(input_file,
                           chrom,
                           quality_threshold,
                           start=minimum,
                           stop=maximum,
                           min_length=min_length,
                           max_length=max_length,
                           verbose=(verbose>=2))

    if (verbose):
        stderr.write("Done reading fragments, preparing for WPS calculation.\n")
    # check if no fragments exist on this interval
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
        scores['contig'] = chrom
    else:

        window_centers = np.arange(start, stop, dtype=np.int64)
        window_starts = np.zeros(stop-start)
        window_stops = np.zeros(stop-start)
        np.rint(window_centers - window_size * 0.5, window_starts)
        np.rint(window_centers + window_size * 0.5 - 1, window_stops)
        # stops are inclusive
        # array to store positions and scores
        scores = np.zeros(
            stop-start,
            dtype=[
                ('contig', 'U16'),
                ('start', 'i8'),
                ('wps', 'i8'),
            ]
        )
        for i in range(stop-start):
            scores[i] = _single_nt_wps(
                chrom,
                window_starts[i],
                window_stops[i],
                window_centers[i],
                frag_ends
            )

    # TODO: consider switch-case statements and determine if they
    # shouldn't be used for backwards compatability
    if (type(output_file) == str):   # check if output specified

        if (verbose):
            stderr.write('Writing to output file.\n')

        if output_file.endswith(".wig.gz"): # zipped wiggle
            with gzip.open(output_file, 'wt') as out:
                # declaration line
                out.write(
                    f'fixedStep\tchrom={chrom}\tstart={start}\t'
                    f'step={1}\tspan={stop-start}\n'
                    )
                for score in scores['wps']:
                    out.write(f'{score}\n')

        elif output_file.endswith(".wig"):  # wiggle
            with open(output_file, 'wt') as out:
                # declaration line
                out.write(
                    f'fixedStep\tchrom={chrom}\tstart={start}\tstep='
                    f'{1}\tspan={stop-start}\n'
                    )
                for score in scores['wps']:
                    out.write(f'{score}\n')

        elif output_file == '-':    #stdout
            stdout.write(
                f'fixedStep\tchrom={chrom}\tstart={start}\tstep='
                f'{1}\tspan={stop-start}\n'
                )
            for score in scores['wps']:
                stdout.write(f'{score}\n')
            stdout.flush()

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
        stderr.write(f'wps took {end_time - start_time} s to complete\n')

    return scores
