"""
5' end-motif features (Zhou et al., 2023): genome-wide, region, and
interval-stratified k-mer end-motif frequencies plus motif diversity score.
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
    FPROFILE_PATH,
    MIN_QUALITY,
    _genome_window_args,
    _MotifFreqs,
    _MotifsIntervals,
    aggregate_genome_motifs,
    aggregate_interval_motifs,
    resolve_motif_aliases,
    write_motif_freqs,
)

__all__ = [
    "EndMotifFreqs",
    "EndMotifsIntervals",
    "region_end_motifs",
    "end_motifs",
    "interval_end_motifs",
    "FPROFILE_PATH",
    "MIN_QUALITY",
]


class EndMotifFreqs(_MotifFreqs):
    """Genome-wide 5' end-motif k-mer frequencies (Zhou et al., 2023)."""


class EndMotifsIntervals(_MotifsIntervals):
    """Interval-stratified 5' end-motif k-mer counts."""


def region_end_motifs(
    input_file: str,
    contig: str,
    start: int,
    stop: int,
    refseq_file: str | Path,
    k: int = 4,
    fraction_low: int | None = 50,
    fraction_high: int | None = None,
    both_strands: bool = True,
    negative_strand: bool = False,
    output_file: str | None = None,
    quality_threshold: int = MIN_QUALITY,
    verbose: bool | int = False,
) -> dict:
    """Count 5' end-motif k-mers for fragments in a region.

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
        End-motif k-mer length (default 4).
    fraction_low, fraction_high : int, optional
        Minimum/maximum fragment length.
    both_strands : bool, optional
        Use 5' ends of both mates; if ``False`` only the forward end is used
        unless ``negative_strand`` is set.
    negative_strand : bool, optional
        With ``both_strands=False``, use only negative-strand 5' ends.
    output_file : str, optional
        Ignored (kept for signature compatibility).
    quality_threshold : int, optional
        Minimum mapping quality (default 20).
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

    # fraction_low < k breaks py2bit sequence queries.
    if fraction_low < k:
        warnings.warn(
            f"fraction_low={fraction_low} < k={k}, which may cause errors. "
            "Automatically setting fraction_low=k."
        )
        fraction_low = k

    with AlignmentWrapper(
        input_file,
        reference_file=refseq_file,
        quality_threshold=quality_threshold,
    ) as wrapper:
        frag_ends = wrapper.fetch(contig, start, stop)

        kmer_list = gen_kmers(k, "ACGT")
        end_motif_counts = dict(zip(kmer_list, 4**k * [0]))

        with ReferenceWrapper(str(refseq_file), use_lock=False) as refseq:
            chroms_dict = refseq.chroms
            if both_strands:
                for frag in frag_ends:
                    try:
                        forward_kmer = refseq.sequence(
                            contig, int(frag.start), int(frag.start + k)
                        )
                        if len(forward_kmer) == k and "N" not in forward_kmer:
                            end_motif_counts[forward_kmer] += 1
                    except ValueError:
                        continue

                    try:
                        reverse_kmer = refseq.sequence(
                            contig, int(frag.stop - k), int(frag.stop)
                        )
                        if len(reverse_kmer) == k and "N" not in reverse_kmer:
                            end_motif_counts[reverse_complement(reverse_kmer)] += 1
                    except ValueError as e:
                        raise RuntimeError(
                            "Error querying sequence at "
                            f"{contig}:{int(frag.stop - k)}-{int(frag.stop)}. "
                            f"Chrom length: {chroms_dict.get(contig, 'unknown')}. "
                            "Please verify that the reference file matches the "
                            f"fragment file. Original error: {e}"
                        )
            else:
                for frag in frag_ends:
                    if frag.is_forward and not negative_strand:
                        try:
                            forward_kmer = refseq.sequence(
                                contig, int(frag.start), int(frag.start + k)
                            )
                            if len(forward_kmer) == k and "N" not in forward_kmer:
                                end_motif_counts[forward_kmer] += 1
                        except ValueError:
                            continue
                    elif negative_strand:
                        try:
                            reverse_kmer = refseq.sequence(
                                contig, int(frag.stop - k), int(frag.stop)
                            )
                            if len(reverse_kmer) == k and "N" not in reverse_kmer:
                                end_motif_counts[
                                    reverse_complement(reverse_kmer)
                                ] += 1
                        except ValueError:
                            if verbose > 1:
                                stderr.write(
                                    f"Attempt to read interval at {contig}:"
                                    f"{int(frag.stop - k)}-{int(frag.stop)} failed."
                                    "Skipping."
                                )
                            continue

    if verbose:
        stop_time = time()
        stderr.write(
            f"region_end_motifs took {stop_time - start_time} seconds to run\n"
        )

    return end_motif_counts


def _region_end_motifs_star(args) -> NDArray[np.float64]:
    return np.array(list(region_end_motifs(*args).values()), dtype="<f8")


def _region_end_motifs_dict_star(args) -> dict:
    return region_end_motifs(*args)


def end_motifs(
    input_file: str,
    refseq_file: str | Path,
    k: int = 4,
    min_length: int | None = 50,
    max_length: int | None = None,
    both_strands: bool = True,
    negative_strand: bool = False,
    output_file: None | str = None,
    quality_threshold: int = 30,
    workers: int = 1,
    verbose: bool | int = False,
    fraction_low: int | None = None,
    fraction_high: int | None = None,
) -> EndMotifFreqs:
    """Compute genome-wide 5' end-motif frequencies (Zhou et al., 2023).

    Parameters
    ----------
    input_file : str
        BAM/CRAM/fragment input.
    refseq_file : str or Path
        Reference genome (``.2bit`` or FASTA).
    k : int, optional
        End-motif k-mer length (default 4).
    min_length, max_length : int, optional
        Fragment-length filter (default min 50).
    both_strands : bool, optional
        Use both mates' 5' ends (default ``True``).
    negative_strand : bool, optional
        With ``both_strands=False``, use only negative-strand ends.
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
    EndMotifFreqs
        Normalized genome-wide end-motif frequencies.
    """
    if verbose:
        start_time = time()

    min_length, max_length = resolve_motif_aliases(
        min_length, max_length, fraction_low, fraction_high
    )

    if min_length is not None and min_length < k:
        warnings.warn(
            f"min_length={min_length} < k={k}, which may cause errors. "
            "Automatically setting min_length=k."
        )
        min_length = k

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
        _region_end_motifs_star,
        EndMotifFreqs,
        k,
        quality_threshold,
        workers,
        verbose,
        "Reading 1mb windows",
        "Counting end-motifs",
    )

    write_motif_freqs(results, output_file)

    if verbose:
        stop_time = time()
        stdout.write(f"end_motifs took {stop_time - start_time} seconds to run\n")

    return results


def interval_end_motifs(
    input_file: str,
    refseq_file: str | Path,
    intervals: str | Iterable[tuple[str, int, int, str]],
    k: int = 4,
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
) -> EndMotifsIntervals:
    """Compute interval-stratified 5' end-motif frequencies.

    Parameters
    ----------
    input_file : str
        BAM/CRAM/fragment input.
    refseq_file : str or Path
        Reference genome.
    intervals : str or list of tuple
        BED path or list of ``(chrom, start, stop, name)`` tuples.
    k : int, optional
        End-motif k-mer length (default 4).
    min_length, max_length : int, optional
        Fragment-length filter.
    both_strands : bool, optional
        Use both mates' 5' ends (default ``True``).
    negative_strand : bool, optional
        With ``both_strands=False``, use only negative-strand ends.
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
    EndMotifsIntervals
        Per-interval end-motif counts.
    """
    if verbose:
        start_time = time()

    min_length, max_length = resolve_motif_aliases(
        min_length, max_length, fraction_low, fraction_high
    )

    if min_length is not None and min_length < k:
        warnings.warn(
            f"min_length={min_length} < k={k}, which may cause errors. "
            "Automatically setting min_length=k."
        )
        min_length = k

    results = aggregate_interval_motifs(
        input_file,
        refseq_file,
        intervals,
        _region_end_motifs_dict_star,
        EndMotifsIntervals,
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
        stdout.write(f"end_motifs took {stop_time - start_time} seconds to run\n")

    return results


def _cli_mds(file_path: str, sep: str = "\t", header: int = 0) -> None:
    """CLI: print the motif diversity score of a k-mer frequency file."""
    # 30 is a placeholder quality threshold; it does not affect MDS.
    motifs = EndMotifFreqs.from_file(file_path, 30, sep, header)
    stdout.write(f"{motifs.motif_diversity_score()}\n")


def _cli_regional_mds(
    file_path: str, file_out: str, sep: str = ",", header: int = 0
) -> None:
    """CLI: write the regional MDS (rMDS) of each region to ``file_out``."""
    motifs = EndMotifsIntervals.from_file(file_path, 30, sep, header)
    motifs.mds_bed(file_out)
