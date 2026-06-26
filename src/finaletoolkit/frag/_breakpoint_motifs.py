"""
Breakpoint-motif features: k-mers centered on fragment start/stop breakpoints,
genome-wide and interval-stratified, plus motif diversity score.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from sys import stderr, stdout
from time import time
from typing import Iterable

import numpy as np
from numpy.typing import NDArray

from finaletoolkit.io.alignment import AlignmentWrapper
from finaletoolkit.io.reference import ReferenceWrapper
from finaletoolkit.utils import gen_kmers, reverse_complement

from ._motif_common import (
    _genome_window_args,
    _MotifFreqs,
    _MotifsIntervals,
    aggregate_genome_motifs,
    aggregate_interval_motifs,
    resolve_motif_aliases,
    write_motif_freqs,
)

__all__ = [
    "BreakpointMotifFreqs",
    "BreakpointMotifsIntervals",
    "region_breakpoint_motifs",
    "breakpoint_motifs",
    "interval_breakpoint_motifs",
]


class BreakpointMotifFreqs(_MotifFreqs):
    """Genome-wide breakpoint-motif k-mer frequencies."""

    def __init__(self, kmer_frequencies, k, quality_threshold: int = 30) -> None:
        super().__init__(kmer_frequencies, k, quality_threshold)


class BreakpointMotifsIntervals(_MotifsIntervals):
    """Interval-stratified breakpoint-motif k-mer counts."""

    def __init__(self, intervals, k, quality_threshold: int = 30) -> None:
        super().__init__(intervals, k, quality_threshold)


def region_breakpoint_motifs(
    input_file: str,
    contig: str,
    start: int,
    stop: int,
    refseq_file: str | Path,
    k: int = 6,
    fraction_low: int = 10,
    fraction_high: int = 600,
    both_strands: bool = True,
    negative_strand: bool = False,
    output_file: str | None = None,
    quality_threshold: int = 30,
    verbose: bool | int = False,
) -> dict:
    """Count breakpoint-motif k-mers for fragments in a region.

    Unlike end motifs, the k-mer is read symmetrically around each breakpoint
    (``k//2`` bp on each side of the fragment start/stop).

    Parameters
    ----------
    input_file : str
        BAM/CRAM/fragment input.
    contig : str
        Contig name.
    start : int
        0-based start coordinate.
    stop : int
        1-based stop coordinate.
    refseq_file : str or Path
        Reference the input was aligned to.
    k : int, optional
        Breakpoint-motif k-mer length (default 6).
    fraction_low, fraction_high : int, optional
        Minimum/maximum fragment length (defaults 10/600).
    both_strands : bool, optional
        Use both breakpoints; if ``False`` only the forward one is used unless
        ``negative_strand`` is set.
    negative_strand : bool, optional
        With ``both_strands=False``, use only the negative-strand breakpoint.
    output_file : str, optional
        Ignored (kept for signature compatibility).
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    verbose : bool or int, optional
        Print timing.

    Returns
    -------
    dict
        Mapping of k-mer to count (all ``4**k`` k-mers present).
    """
    if verbose:
        start_time = time()

    if both_strands and negative_strand:
        raise ValueError("Cannot have both both_strands and negative_strand.")

    half_k = k // 2

    with AlignmentWrapper(
        input_file,
        reference_file=refseq_file,
        quality_threshold=quality_threshold,
    ) as wrapper:
        frag_ends = wrapper.fetch(contig, start, stop)

        kmer_list = gen_kmers(k, "ACGT")
        breakpoint_motif_counts = dict(zip(kmer_list, 4**k * [0]))

        with ReferenceWrapper(str(refseq_file), use_lock=False) as refseq:
            chroms_dict = refseq.chroms
            for frag in frag_ends:
                # Skip fragments whose breakpoint window falls off the contig.
                if (frag.start - half_k) < 0 or (
                    frag.start + half_k
                ) >= chroms_dict.get(frag.contig, 0):
                    warnings.warn(
                        f"Fragment {frag.contig}:{frag.start}-{frag.stop} is "
                        "too close to the end of chromosome. Skipping."
                    )
                    continue

                use_forward = both_strands or (
                    frag.is_forward and not negative_strand
                )
                use_reverse = both_strands or negative_strand

                if use_forward:
                    try:
                        forward_kmer = refseq.sequence(
                            contig, frag.start - half_k, frag.start + half_k
                        )
                        if len(forward_kmer) != k:
                            warnings.warn(
                                f"Skipped fragment {contig}:{frag.start}-"
                                f"{frag.stop} due to length discrepancy with "
                                "motif. This may be due to the fragment being "
                                "aligned to the end of a mitochondrial DNA."
                            )
                            continue
                        if "N" not in forward_kmer:
                            breakpoint_motif_counts[forward_kmer] += 1
                    except ValueError:
                        continue

                if use_reverse:
                    try:
                        reverse_kmer = refseq.sequence(
                            contig, frag.stop - half_k, frag.stop + half_k
                        )
                        if len(reverse_kmer) != k:
                            warnings.warn(
                                f"Skipped fragment {contig}:{frag.start}-"
                                f"{frag.stop} due to length discrepancy with "
                                "motif. This may be due to the fragment being "
                                "aligned to the end of a mitochondrial DNA."
                            )
                            continue
                        if "N" not in reverse_kmer:
                            breakpoint_motif_counts[
                                reverse_complement(reverse_kmer)
                            ] += 1
                    except (RuntimeError, ValueError):
                        if verbose > 1:
                            stderr.write(
                                f"Attempt to read interval at {contig}:"
                                f"{int(frag.stop - half_k)}-"
                                f"{int(frag.stop + half_k)} failed.Skipping."
                            )
                        continue

    if verbose:
        stop_time = time()
        stderr.write(
            f"region_breakpoint_motifs took {stop_time - start_time} seconds "
            "to run\n"
        )

    return breakpoint_motif_counts


def _region_breakpoint_motifs_star(args) -> NDArray[np.float64]:
    return np.array(list(region_breakpoint_motifs(*args).values()), dtype="<f8")


def _region_breakpoint_motifs_dict_star(args) -> dict:
    return region_breakpoint_motifs(*args)


def breakpoint_motifs(
    input_file: str,
    refseq_file: str | Path,
    k: int = 6,
    min_length: int = 50,
    max_length: int = None,
    both_strands: bool = True,
    negative_strand: bool = False,
    output_file: None | str = None,
    quality_threshold: int = 30,
    workers: int = 1,
    verbose: bool | int = False,
    fraction_low: int | None = None,
    fraction_high: int | None = None,
) -> BreakpointMotifFreqs:
    """Compute genome-wide breakpoint-motif frequencies.

    Parameters
    ----------
    input_file : str
        BAM/CRAM/fragment input.
    refseq_file : str or Path
        Reference genome.
    k : int, optional
        Breakpoint-motif k-mer length (default 6).
    min_length, max_length : int, optional
        Fragment-length filter (default min 50).
    both_strands : bool, optional
        Use both breakpoints (default ``True``).
    negative_strand : bool, optional
        With ``both_strands=False``, use only negative-strand breakpoints.
    output_file : str, optional
        TSV/CSV output path.
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    workers : int, optional
        Worker-process count (default 1).
    verbose : bool or int, optional
        Print progress/timing.
    fraction_low, fraction_high : int, optional
        Deprecated aliases for ``min_length``/``max_length``.

    Returns
    -------
    BreakpointMotifFreqs
        Normalized genome-wide breakpoint-motif frequencies.
    """
    if verbose:
        start_time = time()

    min_length, max_length = resolve_motif_aliases(
        min_length, max_length, fraction_low, fraction_high
    )

    with ReferenceWrapper(str(refseq_file), use_lock=False) as refseq:
        chroms = refseq.chroms

    intervals = _genome_window_args(
        input_file,
        refseq_file,
        chroms,
        k,
        min_length,
        max_length,
        both_strands,
        negative_strand,
        quality_threshold,
        verbose,
    )

    results = aggregate_genome_motifs(
        intervals,
        _region_breakpoint_motifs_star,
        BreakpointMotifFreqs,
        k,
        quality_threshold,
        workers,
        verbose,
        "Reading 1mb windows",
        "Counting breakpoint-motifs",
    )

    write_motif_freqs(results, output_file)

    if verbose:
        stop_time = time()
        stdout.write(
            f"breakpoint_motifs took {stop_time - start_time} seconds to run\n"
        )

    return results


def interval_breakpoint_motifs(
    input_file: str,
    refseq_file: str | Path,
    intervals: str | Iterable[tuple[str, int, int, str]],
    k: int = 6,
    min_length: int | None = 50,
    max_length: int | None = None,
    both_strands: bool = True,
    negative_strand: bool = False,
    output_file: str | None = None,
    quality_threshold: int = 30,
    workers: int = 1,
    verbose: bool | int = False,
    fraction_low: int | None = None,
    fraction_high: int | None = None,
) -> BreakpointMotifsIntervals:
    """Compute interval-stratified breakpoint-motif frequencies.

    Parameters
    ----------
    input_file : str
        BAM/CRAM/fragment input.
    refseq_file : str or Path
        Reference genome.
    intervals : str or list of tuple
        BED path or list of ``(chrom, start, stop, name)`` tuples.
    k : int, optional
        Breakpoint-motif k-mer length (default 6).
    min_length, max_length : int, optional
        Fragment-length filter.
    both_strands : bool, optional
        Use both breakpoints (default ``True``).
    negative_strand : bool, optional
        With ``both_strands=False``, use only negative-strand breakpoints.
    output_file : str, optional
        TSV/CSV output path.
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    workers : int, optional
        Worker-process count (default 1).
    verbose : bool or int, optional
        Print progress/timing.
    fraction_low, fraction_high : int, optional
        Deprecated aliases for ``min_length``/``max_length``.

    Returns
    -------
    BreakpointMotifsIntervals
        Per-interval breakpoint-motif counts.
    """
    if verbose:
        start_time = time()

    min_length, max_length = resolve_motif_aliases(
        min_length, max_length, fraction_low, fraction_high
    )

    results = aggregate_interval_motifs(
        input_file,
        refseq_file,
        intervals,
        _region_breakpoint_motifs_dict_star,
        BreakpointMotifsIntervals,
        k,
        min_length,
        max_length,
        both_strands,
        negative_strand,
        quality_threshold,
        workers,
        verbose,
        "Reading intervals",
    )

    write_motif_freqs(results, output_file)

    if verbose:
        stop_time = time()
        stdout.write(
            f"breakpoint_motifs took {stop_time - start_time} seconds to run\n"
        )

    return results


def _cli_mds(file_path: str, sep: str = "\t", header: int = 0) -> None:
    """CLI: print the motif diversity score of a k-mer frequency file."""
    motifs = BreakpointMotifFreqs.from_file(file_path, 30, sep, header)
    stdout.write(f"{motifs.motif_diversity_score()}\n")


def _cli_regional_mds(
    file_path: str, file_out: str, sep: str = ",", header: int = 0
) -> None:
    """CLI: write the regional MDS (rMDS) of each region to ``file_out``."""
    motifs = BreakpointMotifsIntervals.from_file(file_path, 30, sep, header)
    motifs.mds_bed(file_out)
