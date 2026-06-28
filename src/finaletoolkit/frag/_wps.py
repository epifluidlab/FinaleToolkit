"""
Windowed Protection Score (WPS) over a single interval (Snyder et al., 2016).
"""
from __future__ import annotations

import gzip
import time
import warnings
from pathlib import Path
from sys import stderr, stdout
from typing import Union

import numpy as np
import pysam
from numba import jit

from finaletoolkit.utils import frag_array

__all__ = ["wps"]

# Structured dtype for per-position WPS scores.
_WPS_DTYPE = [("contig", "U16"), ("start", "i8"), ("wps", "i8")]


@jit(nopython=True)
def _single_nt_wps(
    chrom: str,
    window_start: int,
    window_stop: int,
    window_position: int,
    frag_ends: np.ndarray,
) -> tuple:
    """Compute the WPS at one position: spanning fragments minus fragment ends.

    A fragment *spans* the window if it starts before and ends after it; a
    fragment *ends in* the window if either endpoint lies within
    ``[window_start, window_stop]``.
    """
    is_spanning = (frag_ends["start"] < window_start) * (
        frag_ends["stop"] > window_stop
    )
    num_spanning = np.sum(is_spanning)

    is_start_in = (frag_ends["start"] >= window_start) * (
        frag_ends["start"] <= window_stop
    )
    is_stop_in = (frag_ends["stop"] >= window_start) * (
        frag_ends["stop"] <= window_stop
    )
    is_end_in = np.logical_or(is_start_in, is_stop_in)
    num_end_in = np.sum(is_end_in)

    return (chrom, window_position, num_spanning - num_end_in)


def wps(
    input_file: Union[str, pysam.AlignmentFile],
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
    reference_file: str | Path | None = None,
) -> np.ndarray:
    """Compute raw Windowed Protection Scores over ``chrom:[start, stop)``.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM/CRAM/fragment input.
    chrom : str
        Contig name.
    start, stop : int
        Region bounds (0-based, half-open).
    chrom_size : int
        Length of ``chrom`` (bounds the fragment-fetch window).
    output_file : str, optional
        ``.wig``/`.wig.gz` path or ``"-"`` for stdout.
    window_size : int, optional
        WPS window width (default 120, i.e. L-WPS).
    min_length, max_length : int, optional
        Fragment-length filter (defaults 120/180, i.e. long fraction).
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    verbose : bool or int, optional
        Print progress/timing.
    fraction_low, fraction_high : int, optional
        Deprecated aliases for ``min_length``/``max_length``.
    reference_file : str or Path, optional
        Reference genome (required for CRAM).

    Returns
    -------
    numpy.ndarray
        Structured array with fields ``('contig', 'start', 'wps')``.  Empty for
        a degenerate interval.
    """
    if verbose:
        start_time = time.time()
        stderr.write("[finaletoolkit-wps] Reading fragments\n")
        stderr.write(f"Region: {chrom}:{start}-{stop}\n")

    # Resolve deprecated aliases (preserves the original conflict behavior:
    # supplying both spellings is an error).
    if fraction_low is not None and min_length is None:
        min_length = fraction_low
        warnings.warn(
            "fraction_low is deprecated. Use min_length instead.",
            category=DeprecationWarning,
            stacklevel=2,
        )
    elif fraction_low is not None and min_length is not None:
        warnings.warn(
            "fraction_low is deprecated. Use min_length instead.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        raise ValueError("fraction_low and min_length cannot both be specified")

    if fraction_high is not None and max_length is None:
        max_length = fraction_high
        warnings.warn(
            "fraction_high is deprecated. Use max_length instead.",
            category=DeprecationWarning,
            stacklevel=2,
        )
    elif fraction_high is not None and max_length is not None:
        warnings.warn(
            "fraction_high is deprecated. Use max_length instead.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        raise ValueError("fraction_high and max_length cannot both be specified")

    start = int(start)
    stop = int(stop)

    if stop <= start:
        warnings.warn(
            f"[wps] {chrom}:{start}-{stop} is a degenerate interval "
            "(stop <= start); skipping.",
            UserWarning,
            stacklevel=2,
        )
        return np.zeros(0, dtype=_WPS_DTYPE)

    # Fetch fragments from a padded window so fragments spanning the edges of
    # [start, stop) are included.
    minimum = max(round(start - max_length), 0)
    maximum = min(round(stop + max_length), chrom_size)

    frag_ends = frag_array(
        input_file,
        chrom,
        quality_threshold,
        start=minimum,
        stop=maximum,
        min_length=min_length,
        max_length=max_length,
        verbose=(verbose >= 2),
        reference_file=reference_file,
    )

    if verbose:
        stderr.write("Done reading fragments, preparing for WPS calculation.\n")

    # The per-position loop below also yields the correct all-zero result when
    # no fragments are present, so no separate empty-interval branch is needed.
    window_centers = np.arange(start, stop, dtype=np.int64)
    window_starts = np.rint(window_centers - window_size * 0.5)
    window_stops = np.rint(window_centers + window_size * 0.5 - 1)  # inclusive

    scores = np.zeros(stop - start, dtype=_WPS_DTYPE)
    for i in range(stop - start):
        scores[i] = _single_nt_wps(
            chrom,
            window_starts[i],
            window_stops[i],
            window_centers[i],
            frag_ends,
        )

    if isinstance(output_file, str):
        if verbose:
            stderr.write("Writing to output file.\n")
        _write_wig(output_file, chrom, start, stop, scores)
    elif output_file is not None:
        raise TypeError(
            f'output_file is unsupported type "{type(input_file)}". '
            "output_file should be a string specifying the path of the file "
            "to output scores to."
        )

    if verbose:
        end_time = time.time()
        stderr.write(f"wps took {end_time - start_time} s to complete\n")

    return scores


def _write_wig(output_file, chrom, start, stop, scores) -> None:
    """Write per-position WPS scores in fixedStep WIG format."""
    header = (
        f"fixedStep\tchrom={chrom}\tstart={start}\tstep={1}\tspan={stop - start}\n"
    )
    if output_file.endswith(".wig.gz"):
        with gzip.open(output_file, "wt") as out:
            out.write(header)
            for score in scores["wps"]:
                out.write(f"{score}\n")
    elif output_file.endswith(".wig"):
        with open(output_file, "wt") as out:
            out.write(header)
            for score in scores["wps"]:
                out.write(f"{score}\n")
    elif output_file == "-":
        stdout.write(header)
        for score in scores["wps"]:
            stdout.write(f"{score}\n")
        stdout.flush()
    else:
        raise ValueError("output_file can only have suffixes .wig or .wig.gz.")
