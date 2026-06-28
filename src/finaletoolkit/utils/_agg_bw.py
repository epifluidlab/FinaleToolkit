"""
Aggregate a bigWig signal across strand-oriented intervals.
"""
from __future__ import annotations

import gzip
import time
from os import PathLike
from sys import stderr
from typing import Union

import numpy as np
import pyBigWig as pbw

__all__ = ["agg_bw"]


def agg_bw(
    input_file: Union[str, PathLike],
    interval_file: Union[str, PathLike],
    output_file: Union[str, PathLike],
    median_window_size: int = 1,
    mean: bool = False,
    verbose: bool = False,
) -> np.ndarray:
    """Aggregate bigWig signal over BED intervals (strand-aware).

    The median filter used by ``adjust-wps`` trims each interval by half the
    window size.  Account for that either by supplying smaller intervals or by
    passing the original ``median_window_size`` here (do not do both).

    Parameters
    ----------
    input_file : str or path
        bigWig of per-base signal.
    interval_file : str or path
        BED of intervals; column 6 must contain the strand.
    output_file : str or path
        Output WIG path.
    median_window_size : int, optional
        Filter window used upstream (default 1 = no trimming; 120 replicates
        Snyder et al.).
    mean : bool, optional
        Divide the aggregate by the number of intervals (mean instead of sum).
    verbose : bool, optional
        Print progress/timing.

    Returns
    -------
    numpy.ndarray
        The aggregated per-position signal.

    Raises
    ------
    ValueError
        If ``interval_file`` is not BED or ``output_file`` is not ``.wig``.
    """
    if verbose:
        start_time = time.time()
        stderr.write("Reading intervals from bed...\n")

    if not (
        str(interval_file).endswith(".bed")
        or str(interval_file).endswith(".bed.gz")
    ):
        raise ValueError("Invalid filetype for interval_file.")

    intervals = []
    opener = gzip.open if str(interval_file).endswith(".gz") else open
    with opener(interval_file, "rt") as file:
        for line in file:
            contents = line.split("\t")
            contig = contents[0]
            start = int(contents[1])
            stop = int(contents[2])
            strand = contents[5]
            intervals.append((contig, int(start), int(stop), strand.strip()))

    with pbw.open(str(input_file), "r") as raw_wps:  # pbw needs str, not Path
        # Interval length after trimming by the median-filter window.
        interval_size = intervals[0][2] - intervals[0][1] - median_window_size
        agg_scores = np.zeros(interval_size, dtype=np.int64)
        num_intervals_added = 0
        for contig, start, stop, strand in intervals:
            try:
                signal = raw_wps.values(contig, start, stop)
                if signal is None:
                    print(
                        "There was no information found in the interval: ",
                        contig,
                        start,
                        stop,
                    )
                    continue
                values = np.nan_to_num(np.array(signal), nan=0)
            except RuntimeError as e:
                print(e)
                continue

            # Trim the ends removed by the upstream median filter.
            trimmed = values[median_window_size // 2 : -median_window_size // 2]
            if trimmed.shape[0] != interval_size:
                print(
                    f"Trimmed size {trimmed.shape[0]} for {contig}:{start}"
                    f"-{stop} is not equal to "
                    f"interval size {interval_size}. Skipping."
                )
                continue

            if strand == "+":
                agg_scores = agg_scores + trimmed
                num_intervals_added += 1
            elif strand == "-":
                agg_scores = agg_scores + np.flip(trimmed)
                num_intervals_added += 1
            elif verbose:
                stderr.write(
                    "A segment without strand was encountered. Skipping."
                )

    if mean:
        agg_scores = agg_scores / num_intervals_added

    if str(output_file).endswith("wig"):
        with open(output_file, "wt") as out:
            if verbose:
                stderr.write("File opened! Writing...\n")
            out.write(
                f"fixedStep\tchrom=.\tstart={-interval_size // 2}\tstep={1}\t"
                f"span={interval_size}\n"
            )
            for score in agg_scores:
                out.write(f"{score}\n")
    else:
        raise ValueError(
            "The output_file is an unaccepted type. Must be a wiggle file "
            "ending in .wig"
        )

    if verbose:
        end_time = time.time()
        stderr.write(
            f"Aggregating bigWig took {end_time - start_time} s to run.\n"
        )

    return agg_scores
