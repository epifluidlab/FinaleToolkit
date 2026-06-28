"""
Fragment coverage over intervals.

Coverage is estimated by counting fragments (by the configured intersect
policy) that fall in each interval, optionally scaled and/or normalized by the
genome-wide fragment count.
"""
from __future__ import annotations

import gzip
import sys
import time
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import NamedTuple, Union

import pysam
from tqdm import tqdm

from finaletoolkit.utils import frag_generator, get_intervals

__all__ = ["coverage", "single_coverage", "CoverageResult"]


class CoverageResult(NamedTuple):
    """Coverage over a single interval.

    A named tuple: it unpacks and indexes like
    ``(contig, start, stop, name, coverage)`` and also exposes named fields.

    Attributes
    ----------
    contig : str or None
        Interval contig.
    start : int or None
        0-based start coordinate.
    stop : int or None
        Stop coordinate.
    name : str
        Interval name.
    coverage : float
        Coverage value over the interval.
    """

    contig: str | None
    start: int | None
    stop: int | None
    name: str
    coverage: float


def single_coverage(
    input_file: Union[str, pysam.TabixFile, pysam.AlignmentFile, Path],
    contig: str | None = None,
    start: int | None = 0,
    stop: int | None = None,
    name: str | None = ".",
    min_length: int | None = None,
    max_length: int | None = None,
    intersect_policy: str = "midpoint",
    quality_threshold: int = 30,
    verbose: bool | int = False,
    reference_file: str | Path | None = None,
) -> CoverageResult:
    """Estimate fragment coverage over a single region.

    Counts fragments (per ``intersect_policy``) falling in
    ``contig:[start, stop)``.  Not suitable when fragment sizes approach the
    region size.

    Parameters
    ----------
    input_file : str or pysam handle
        BAM/CRAM/fragment input.
    contig : str, optional
        Contig name.
    start : int, optional
        0-based start (default 0).
    stop : int, optional
        1-based stop (default: end of contig).
    name : str, optional
        Interval name (``'.'`` if ``None``).
    min_length, max_length : int, optional
        Fragment-length filter.
    intersect_policy : {"midpoint", "any"}, optional
        Region-membership policy (default ``"midpoint"``).
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    verbose : bool or int, optional
        Print timing information.
    reference_file : str or Path, optional
        Reference genome (required for CRAM).

    Returns
    -------
    CoverageResult
        ``(contig, start, stop, name, coverage)``.
    """
    if verbose:
        start_time = time.time()
        tqdm.write(
            f"""
            input_file: {input_file}
            contig: {contig}
            start: {start}
            stop: {stop}
            name: {name}
            min_length: {min_length}
            max_length: {max_length}
            intersect_policy: {intersect_policy}
            quality_threshold: {quality_threshold}
            verbose: {verbose}
            """
        )

    coverage = 0
    frags = frag_generator(
        input_file=input_file,
        contig=contig,
        quality_threshold=quality_threshold,
        start=start,
        stop=stop,
        min_length=min_length,
        max_length=max_length,
        intersect_policy=intersect_policy,
        reference_file=reference_file,
    )
    for _ in frags:
        coverage += 1

    if verbose:
        end_time = time.time()
        tqdm.write(f"single_coverage took {end_time - start_time} s to complete\n")

    adjusted_name = "." if name is None else name
    return CoverageResult(contig, start, stop, adjusted_name, coverage)


def _single_coverage_star(partial_coverage, interval) -> CoverageResult:
    contig, start, stop, name = interval
    return partial_coverage(contig=contig, start=start, stop=stop, name=name)


def coverage(
    input_file: Union[str, pysam.TabixFile, pysam.AlignmentFile, Path],
    interval_file: str,
    output_file: str,
    scale_factor: float = 1.0,
    min_length: int | None = None,
    max_length: int | None = None,
    normalize: bool = False,
    intersect_policy: str = "midpoint",
    quality_threshold: int = 30,
    workers: int = 1,
    verbose: Union[bool, int] = False,
    reference_file: str | Path | None = None,
) -> list[CoverageResult]:
    """Estimate fragment coverage over every interval in a BED file.

    Parameters
    ----------
    input_file : str or pysam handle
        BAM/CRAM/fragment input.
    interval_file : str
        BED4 file of intervals.
    output_file : str
        Output BED/`.bedgraph`/`.bed.gz` path, or ``"-"`` for stdout.  ``None``
        suppresses file output and only returns results.
    scale_factor : float, optional
        Multiplier applied to coverage values (default 1.0).
    min_length, max_length : int, optional
        Fragment-length filter.
    normalize : bool, optional
        Divide ``scale_factor`` by the genome-wide coverage before scaling.
    intersect_policy : {"midpoint", "any"}, optional
        Region-membership policy (default ``"midpoint"``).
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    workers : int, optional
        Worker-process count (default 1).
    verbose : bool or int, optional
        Print timing/config information.
    reference_file : str or Path, optional
        Reference genome (required for CRAM).

    Returns
    -------
    list of CoverageResult
        Scaled coverage for each interval.
    """
    return_val: list[CoverageResult] = []
    if verbose:
        start_time = time.time()
        tqdm.write(
            f"""
            input_file: {input_file}
            interval file: {interval_file}
            output_file: {output_file}
            scale_factor: {scale_factor}
            min_length: {min_length}
            max_length: {max_length}
            intersect_policy: {intersect_policy}
            quality_threshold: {quality_threshold}
            workers: {workers}
            normalize: {normalize}
            verbose: {verbose} \n
            """
        )
        tqdm.write("Creating process pool\n")

    pool = Pool(processes=workers)
    try:
        # Kick off the genome-wide total coverage asynchronously if normalizing.
        if normalize:
            total_coverage_results = pool.apply_async(
                single_coverage,
                (input_file, None, 0, None, "."),
                {
                    "min_length": min_length,
                    "max_length": max_length,
                    "intersect_policy": intersect_policy,
                    "quality_threshold": quality_threshold,
                    "verbose": verbose,
                    "reference_file": reference_file,
                },
            )

        intervals = get_intervals(interval_file)

        if verbose:
            tqdm.write("calculating coverage\n")

        partial_single_coverage = partial(
            single_coverage,
            input_file=input_file,
            min_length=min_length,
            max_length=max_length,
            intersect_policy=intersect_policy,
            quality_threshold=quality_threshold,
            verbose=max(0, verbose - 1),
            reference_file=reference_file,
        )
        coverages = pool.imap(
            partial(_single_coverage_star, partial_single_coverage),
            intervals,
            chunksize=max(len(intervals) // workers, 1),
        )

        if normalize:
            total_coverage = total_coverage_results.get()
            if verbose:
                tqdm.write(f"Total coverage is {total_coverage}\n")
            scale_factor /= total_coverage[4]

        output_is_file = False
        if output_file is not None:
            if verbose:
                tqdm.write("Writing results to output\n")
            try:
                if output_file.endswith(".bed") or output_file.endswith(".bedgraph"):
                    output_is_file = True
                    output = open(output_file, "w")
                elif output_file.endswith(".bed.gz"):
                    output = gzip.open(output_file, "wt")
                    output_is_file = True
                elif output_file == "-":
                    output = sys.stdout
                else:
                    raise ValueError(
                        "output_file should have .bed or .bed.gz as suffix"
                    )

                if output_file.endswith(".bedgraph"):
                    for contig, start, stop, name, cov in coverages:
                        output.write(
                            f"{contig}\t{start}\t{stop}\t{cov * scale_factor}\n"
                        )
                        return_val.append(
                            CoverageResult(contig, start, stop, name, cov * scale_factor)
                        )
                else:
                    for contig, start, stop, name, cov in coverages:
                        output.write(
                            f"{contig}\t{start}\t{stop}\t{name}\t"
                            f"{cov * scale_factor}\n"
                        )
                        return_val.append(
                            CoverageResult(contig, start, stop, name, cov * scale_factor)
                        )
            finally:
                if output_is_file:
                    output.close()
        else:
            return_val = [
                CoverageResult(contig, start, stop, name, cov * scale_factor)
                for (contig, start, stop, name, cov) in coverages
            ]
    finally:
        pool.close()

    if verbose:
        end_time = time.time()
        tqdm.write(f"coverage took {end_time - start_time} s to complete\n")
    return return_val
