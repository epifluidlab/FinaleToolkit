"""Cleavage Profiler.

Computes the cleavage profile of Zhou et al. (2022,
https://doi.org/10.1073/pnas.2209852119): the proportion of fragment ends at
each position over the depth at that position (as a percentage).
"""
from __future__ import annotations

import gzip
import time
import warnings
from multiprocessing import Pool
from pathlib import Path
from sys import stderr, stdin
from typing import Union

import numpy as np
import pyBigWig as pbw

from finaletoolkit.utils import (
    chrom_sizes_to_dict,
    chrom_sizes_to_list,
    frag_array,
)
from finaletoolkit.utils.typing import FragFile

__all__ = ["cleavage_profile", "multi_cleavage_profile"]

# Structured dtype for per-position cleavage proportions.
_CLEAVAGE_DTYPE = [("contig", "U16"), ("pos", "i8"), ("proportion", "f8")]


def _coverage_and_ends(
    frag_starts: np.ndarray,
    frag_stops: np.ndarray,
    frag_strands: np.ndarray,
    adj_start: int,
    adj_stop: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-position fragment depth and fragment-end counts over an interval.

    Equivalent to broadcasting each fragment against each position (an
    ``(n_fragments, n_positions)`` boolean matrix) but computed with a
    coverage difference array instead, so memory is ``O(n_positions)``
    rather than ``O(n_fragments * n_positions)``.

    Parameters
    ----------
    frag_starts, frag_stops : numpy.ndarray
        Fragment start/stop coordinates (half-open, ``[start, stop)``).
    frag_strands : numpy.ndarray
        Boolean array, ``True`` for ``+`` strand fragments.
    adj_start, adj_stop : int
        Half-open interval of positions to compute over.

    Returns
    -------
    depth : numpy.ndarray
        Number of fragments covering each position, shape
        ``(adj_stop - adj_start,)``.
    ends : numpy.ndarray
        Number of fragment ends (start of + strand fragments, stop of -
        strand fragments) at each position, same shape as ``depth``.
    """
    n = adj_stop - adj_start

    raw_start_idx = frag_starts - adj_start
    raw_stop_idx = frag_stops - adj_start

    # Depth via a coverage difference array: +1 where a fragment starts,
    # -1 where it stops, then a cumulative sum. Indices are clipped to
    # [0, n] since fragments may extend beyond the interval (frag_array is
    # called with intersect_policy="any").
    depth_diff = np.zeros(n + 1, dtype=np.int64)
    np.add.at(depth_diff, np.clip(raw_start_idx, 0, n), 1)
    np.add.at(depth_diff, np.clip(raw_stop_idx, 0, n), -1)
    depth = np.cumsum(depth_diff[:-1])

    # Fragment ends: forward fragments end at their start (+ strand),
    # reverse fragments end at their stop (- strand). Unlike depth, an
    # end that falls outside the interval must be dropped rather than
    # clipped to the boundary.
    forward_mask = frag_strands
    fwd_idx = raw_start_idx[forward_mask]
    fwd_idx = fwd_idx[(fwd_idx >= 0) & (fwd_idx < n)]
    rev_idx = raw_stop_idx[~forward_mask]
    rev_idx = rev_idx[(rev_idx >= 0) & (rev_idx < n)]
    ends = np.bincount(fwd_idx, minlength=n) + np.bincount(rev_idx, minlength=n)

    return depth, ends


def cleavage_profile(
    input_file: FragFile,
    chrom_size: int,
    contig: str,
    start: int,
    stop: int,
    left: int = 0,
    right: int = 0,
    min_length: int | None = None,
    max_length: int | None = None,
    quality_threshold: int = 30,
    verbose: Union[bool, int] = 0,
    fraction_low: int | None = None,
    fraction_high: int | None = None,
    reference_file: str | Path | None = None,
) -> np.ndarray:
    """Compute the cleavage profile over a single interval.

    Parameters
    ----------
    input_file : str or pysam handle
        BAM/CRAM/fragment input.
    chrom_size : int
        Length of ``contig``.
    contig : str
        Contig name.
    start : int
        0-based start coordinate.
    stop : int
        1-based stop coordinate.
    left, right : int, optional
        Amounts to expand the interval on each side (useful when only CpG
        coordinates are given).
    min_length, max_length : int, optional
        Fragment-length filter.
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
        Structured array with fields ``('contig', 'pos', 'proportion')`` where
        ``proportion`` is a percentage in ``[0, 100]``.
    """
    if verbose:
        start_time = time.time()
        stderr.write(
            f"""
            Calculating cleavage profile
            input_file: {input_file}
            contig: {contig}
            start: {start}
            stop: {stop}
            fraction_low: {fraction_low}
            fraction_high: {fraction_high}
            quality_threshold: {quality_threshold}
            verbose: {verbose}
            """
        )

    # Resolve deprecated aliases (both spellings together is an error).
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

    adj_start = max(start - left, 0)
    adj_stop = min(stop + right, chrom_size)

    frags = frag_array(
        input_file=input_file,
        contig=contig,
        quality_threshold=quality_threshold,
        start=adj_start,
        stop=adj_stop,
        min_length=min_length,
        max_length=max_length,
        intersect_policy="any",
        reference_file=reference_file,
    )

    positions = np.arange(adj_start, adj_stop)
    depth, ends = _coverage_and_ends(
        frags["start"], frags["stop"], frags["strand"], adj_start, adj_stop
    )

    proportions = np.zeros_like(depth, dtype=np.float64)
    non_zero_mask = depth != 0
    proportions[non_zero_mask] = ends[non_zero_mask] / depth[non_zero_mask] * 100

    results = np.zeros_like(proportions, dtype=_CLEAVAGE_DTYPE)
    results["contig"] = contig
    results["pos"] = positions
    results["proportion"] = proportions

    if verbose:
        stderr.write(
            f"cleavage_profile took {time.time() - start_time} s to complete\n"
        )

    return results


def _cleavage_profile_star(args):
    return cleavage_profile(*args)


def multi_cleavage_profile(
    input_file: FragFile,
    interval_file: Union[str, Path],
    chrom_sizes: Union[str, Path],
    left: int = 0,
    right: int = 0,
    min_length: int | None = None,
    max_length: int | None = None,
    quality_threshold: int = 30,
    output_file: str = "-",
    workers: int = 1,
    verbose: Union[bool, int] = 0,
    fraction_low: int | None = None,
    fraction_high: int | None = None,
    reference_file: str | Path | None = None,
) -> str:
    """Compute cleavage profiles over intervals in a (sorted) BED file.

    Parameters
    ----------
    input_file : str or path
        BAM/CRAM/fragment input.
    interval_file : str or path
        Sorted BED of intervals (``"-"`` reads stdin).
    chrom_sizes : str or path
        ``.chrom.sizes`` file (required).
    left, right : int, optional
        Interval expansion applied before merging overlaps.
    min_length, max_length : int, optional
        Fragment-length filter.
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    output_file : str, optional
        ``.bw`` or ``.bed.gz``/`.bedgraph.gz` path (default ``"-"``).
    workers : int, optional
        Worker-process count (default 1).
    verbose : bool or int, optional
        Print progress/timing.
    fraction_low, fraction_high : int, optional
        Deprecated aliases for ``min_length``/``max_length``.
    reference_file : str or Path, optional
        Reference genome (required for CRAM).

    Returns
    -------
    str
        The output path.
    """
    if verbose:
        start_time = time.time()
        stderr.write(
            f"""
            Calculating cleavage profile
            input_file: {input_file}
            interval_file: {interval_file}
            chrom_sizes: {chrom_sizes}
            left: {left}
            right: {right}
            min_length: {min_length}
            max_length: {max_length}
            quality_threshold: {quality_threshold}
            output_file: {output_file}
            workers: {workers}
            verbose: {verbose}
            """
        )

    # Resolve deprecated aliases (both spellings together is an error).
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

    if input_file == "-" and interval_file == "-":
        raise ValueError("input_file and site_bed cannot both read from stdin")

    if chrom_sizes is None:
        raise ValueError("chrom_sizes must be specified.")

    header = chrom_sizes_to_list(chrom_sizes)
    chrom_dict = chrom_sizes_to_dict(chrom_sizes)

    if verbose > 1:
        stderr.write(f"chrom sizes {header}\n")

    contigs, starts, stops = _read_intervals(
        interval_file, left, right, chrom_dict
    )

    size_dict = dict(header)
    sizes = [size_dict[contig] for contig in contigs]
    count = len(contigs)

    if verbose:
        stderr.write("Zipping inputs\n")

    interval_list = zip(
        count * [input_file],
        sizes,
        contigs,
        starts,
        stops,
        count * [0],  # left/right precomputed to avoid double-padding
        count * [0],
        count * [min_length],
        count * [max_length],
        count * [quality_threshold],
        count * [max(verbose - 1, 0)],
        count * [fraction_low],
        count * [fraction_high],
        count * [reference_file],
    )

    if verbose:
        stderr.write("Calculating cleavage profile...\n")

    pool = Pool(workers, maxtasksperchild=500)
    try:
        interval_scores = pool.imap(
            _cleavage_profile_star, interval_list, chunksize=100
        )

        if isinstance(output_file, str):
            if verbose:
                stderr.write(f"Output file {output_file} specified. Opening...\n")
            if output_file.endswith(".bw"):
                _write_bigwig(output_file, header, interval_scores)
            elif (
                output_file.endswith(".bed.gz")
                or output_file.endswith("bedgraph.gz")
                or output_file == "-"
            ):
                _write_bedgraph_gz(output_file, interval_scores)
            else:
                raise ValueError(
                    "output_file can only have suffix .bw, .bedgraph.gz, or "
                    ".bed.gz."
                )
        elif output_file is not None:
            raise TypeError(
                f'output_file is unsupported type "{type(input_file)}". '
                "output_file should be a string specifying the path of the "
                "file to output scores to."
            )
    finally:
        pool.close()

    if verbose:
        end_time = time.time()
        stderr.write(
            f"cleavage profile took {end_time - start_time} s to complete\n"
        )
    return output_file


def _read_intervals(interval_file, left, right, chrom_dict):
    """Parse a sorted BED into merged, expanded intervals."""
    contigs: list[str] = []
    starts: list[int] = []
    stops: list[int] = []

    bed = stdin if interval_file == "-" else open(interval_file)
    try:
        prev_contig = None
        prev_start = 0
        prev_stop = 0
        for line in bed:
            contents = line.split()
            contig = contents[0].strip()
            start, stop = int(contents[1]), int(contents[2])
            if contig not in chrom_dict:
                warnings.warn(
                    f"Skipping interval {contig}:{start}-{stop} from "
                    f"interval_file ({contig} not in chrom_sizes)",
                    UserWarning,
                )
                continue
            start = max(0, start - left)
            stop = min(stop + right, chrom_dict[contig])

            if prev_contig == contig and start < prev_stop:
                prev_stop = max(prev_stop, stop)
            else:
                contigs.append(prev_contig)
                starts.append(prev_start)
                stops.append(prev_stop)
                prev_contig, prev_start, prev_stop = contig, start, stop
        contigs.append(prev_contig)
        starts.append(prev_start)
        stops.append(prev_stop)
    finally:
        if interval_file != "-":
            bed.close()

    # Drop the initial placeholder entry from the prev_* priming.
    return contigs[1:], starts[1:], stops[1:]


def _write_bigwig(output_file, header, interval_scores) -> None:
    """Write per-position cleavage proportions to a bigWig file."""
    with pbw.open(output_file, "w") as bigwig:
        bigwig.addHeader(header)
        last = "None"
        for interval_score in interval_scores:
            contigs = interval_score["contig"]
            starts = interval_score["pos"]
            scores = interval_score["proportion"]

            if contigs.shape == (0,):
                continue
            try:
                bigwig.addEntries(
                    contigs[0],
                    starts[0],
                    values=scores.astype(np.float64),
                    step=1,
                    span=1,
                )
            except RuntimeError as e:
                stderr.write(f"{contigs[0]}:{starts[0]}-{starts[-1] + 1}\n")
                stderr.write(
                    "invalid or out of order interval encountered. "
                    "Skipping to next.\n"
                )
                stderr.write(f"captured error:\n{e}\n")
                stderr.write(f"current output:\n{interval_score}\n")
                stderr.write(f"last output:\n{last}\n")
                continue
            last = interval_score


def _write_bedgraph_gz(output_file, interval_scores) -> None:
    """Write per-position cleavage proportions to a gzip-compressed bedGraph."""
    with gzip.open(output_file, "wt") as bedgraph:
        for interval_score in interval_scores:
            contigs = interval_score["contig"]
            starts = interval_score["pos"]
            scores = interval_score["proportion"]
            stops = starts + 1

            lines = "".join(
                f"{contig}\t{start}\t{stop}\t{score}\n"
                for contig, start, stop, score in zip(contigs, starts, stops, scores)
            )
            bedgraph.write(lines)
