"""
Core fragment/genomic utilities: chrom-size parsing, fragment arrays, quality
filtering, interval parsing, vectorized overlaps, and k-mer helpers.
"""
from __future__ import annotations

import gzip
import itertools
import time
from pathlib import Path
from sys import stderr

import numpy as np
import pysam
from numba import jit
from numpy.typing import NDArray

from ..io.alignment import AlignmentWrapper
from .logging import get_logger
from .typing import ChromSizes, FragFile, Intervals
from ._frag_generator import frag_generator
from ._intervals import (
    _convert_to_list,
    _merge_all_intervals,
    _merge_overlapping_intervals,
    _reduce_overlaps_in_file,
)

logger = get_logger(__name__)

__all__ = [
    "chrom_sizes_to_list",
    "chrom_sizes_to_dict",
    "frag_bam_to_bed",
    "frags_in_region",
    "frag_array",
    "low_quality_read_pairs",
    "_not_read1_or_low_quality",
    "get_intervals",
    "overlaps",
    "gen_kmers",
    "reverse_complement",
    "_merge_overlapping_intervals",
    "_reduce_overlaps_in_file",
    "_convert_to_list",
    "_merge_all_intervals",
]


# -- chrom.sizes parsing ----------------------------------------------------


def chrom_sizes_to_list(chrom_sizes_file: ChromSizes) -> list[tuple[str, int]]:
    """Read a ``.chrom.sizes`` file into a list of ``(name, size)`` tuples.

    Parameters
    ----------
    chrom_sizes_file : str or Path
        Tab-delimited file with contig name and size columns.

    Returns
    -------
    list of (str, int)
        Contig names and sizes, in file order.
    """
    chrom_sizes: list[tuple[str, int]] = []
    with open(chrom_sizes_file, "r") as file:
        for line in file:
            if line != "\n":
                chrom, size = line.strip().split("\t")
                chrom_sizes.append((chrom, int(size)))
    return chrom_sizes


def chrom_sizes_to_dict(chrom_sizes_file: ChromSizes) -> dict[str, int]:
    """Read a ``.chrom.sizes`` file into a ``{name: size}`` dict.

    Parameters
    ----------
    chrom_sizes_file : str or Path
        Tab-delimited file with contig name and size columns.

    Returns
    -------
    dict
        Maps contig name to integer size.
    """
    chrom_sizes: dict[str, int] = {}
    with open(chrom_sizes_file, "r") as file:
        for line in file:
            if line != "\n":
                chrom, size = line.strip().split("\t")
                chrom_sizes[chrom] = int(size)
    return chrom_sizes


# -- fragment export --------------------------------------------------------


def frag_bam_to_bed(
    input_file: str | pysam.AlignmentFile,
    output_file: str,
    contig: str | None = None,
    quality_threshold: int = 30,
    verbose: bool = False,
    reference_file: str | Path | None = None,
) -> None:
    """Write paired-end fragments from a BAM/CRAM file to a 3-column BED.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        Input alignment file.
    output_file : str
        Output BED path (gzip-compressed when it ends with ``.gz``).
    contig : str, optional
        Restrict export to a single contig.
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    verbose : bool, optional
        Print timing information.
    reference_file : str or Path, optional
        Reference genome (required for CRAM).
    """
    if verbose:
        start_time = time.time()
        print("Opening file")

    out = None
    try:
        if output_file.endswith(".gz"):
            out = gzip.open(output_file, "wt")
        else:
            out = open(output_file, "w")

        with AlignmentWrapper(
            input_file,
            reference_file=reference_file,
            quality_threshold=quality_threshold,
        ) as wrapper:
            for frag in wrapper.fetch(contig=contig):
                out.write(f"{frag.contig}\t{frag.start}\t{frag.stop}\n")
    except Exception as e:
        logger.error("An error occurred during BAM to BED conversion: %s", str(e))
    finally:
        if out is not None:
            out.close()

    if verbose:
        end_time = time.time()
        print(
            f"frag_bam_to_bed took {end_time - start_time} s to complete",
            flush=True,
        )


# -- fragment arrays --------------------------------------------------------


@jit(nopython=True)
def frags_in_region(frag_array: NDArray, start: int, stop: int) -> NDArray:
    """Return fragments overlapping ``[start, stop)`` (numba-accelerated).

    A fragment is kept when ``frag.start < stop`` and ``frag.stop >= start``.

    Parameters
    ----------
    frag_array : ndarray
        Structured array with ``'start'`` and ``'stop'`` fields.
    start : int
        Left-most coordinate of the region.
    stop : int
        Right-most coordinate of the region.

    Returns
    -------
    ndarray
        The subset of ``frag_array`` overlapping the region.
    """
    starts = frag_array["start"]
    stops = frag_array["stop"]
    in_region = np.logical_and(np.less(starts, stop), np.greater_equal(stops, start))
    return frag_array[in_region]


def frag_array(
    input_file: FragFile,
    contig: str,
    quality_threshold: int = 30,
    start: int | None = None,
    stop: int | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    intersect_policy: str = "midpoint",
    verbose: bool = False,
    reference_file: str | Path | None = None,
) -> NDArray:
    """Read fragments into a structured ``(start, stop, strand)`` array.

    Parameters
    ----------
    input_file : str or pysam handle
        BAM/CRAM/fragment input.
    contig : str
        Contig to read.
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    start, stop : int, optional
        Region bounds.
    min_length, max_length : int, optional
        Inclusive fragment-length bounds.
    intersect_policy : {"midpoint", "any"}, optional
        Region-membership policy (default ``"midpoint"``).
    verbose : bool, optional
        Print progress to stderr.
    reference_file : str or Path, optional
        Reference genome (required for CRAM).

    Returns
    -------
    ndarray
        1-D structured array with dtype
        ``[('start', 'i8'), ('stop', 'i8'), ('strand', '?')]``.  ``strand`` is
        ``True`` on the ``+`` strand.  Shape ``(0,)`` if no fragments match.
    """
    if verbose:
        stderr.write("frag_array: fetching fragments\n")

    frag_list = [
        (frag_start, frag_stop, strand)
        for _, frag_start, frag_stop, _, strand in frag_generator(
            input_file,
            contig,
            quality_threshold,
            start,
            stop,
            min_length,
            max_length,
            intersect_policy,
            verbose,
            reference_file=reference_file,
        )
    ]

    if verbose:
        stderr.write("frag_array: converting to array\n")

    frag_ends = np.array(
        frag_list, dtype=[("start", "i8"), ("stop", "i8"), ("strand", "?")]
    )

    assert frag_ends.ndim == 1, (
        f"frag_ends has dims {frag_ends.ndim} and shape {frag_ends.shape}"
    )
    return frag_ends


# -- read-level quality filtering -------------------------------------------


def low_quality_read_pairs(read: pysam.AlignedSegment, min_mapq: int = 30) -> bool:
    """Return ``True`` if ``read`` is not a clean, properly-paired alignment.

    Equivalent to samtools ``-F 3852 -f 3`` plus a ``-G 48`` style same-strand
    check and an optional mate mapping-quality (``MQ`` tag) check.  Adapted from
    cofragr's ``frag_summary_in_intervals.py``.

    Parameters
    ----------
    read : pysam.AlignedSegment
        The read to inspect.
    min_mapq : int, optional
        Minimum Phred mapping quality (default 30).

    Returns
    -------
    bool
        ``True`` if the read is low quality, unmapped, or improperly paired.
    """
    if (
        read.is_unmapped  # 0x4
        or read.is_secondary  # 0x100
        or (not read.is_paired)  # not 0x1
        or read.mate_is_unmapped  # 0x8
        or read.is_duplicate  # 0x400
        or read.mapping_quality < min_mapq
        or read.is_qcfail  # 0x200
        or read.is_supplementary  # 0x800
        or (not read.is_proper_pair)  # not 0x2
        or (read.is_reverse and read.mate_is_reverse)  # -G 48
    ):
        return True
    try:
        if read.has_tag("MQ") and read.get_tag("MQ") < min_mapq:
            return True
    except Exception:
        pass

    return False


def _not_read1_or_low_quality(read: pysam.AlignedSegment, min_mapq: int = 30) -> bool:
    """Return ``True`` if ``read`` is not read1 or fails the quality filter."""
    return low_quality_read_pairs(read, min_mapq=min_mapq) or not read.is_read1


# -- interval / overlap helpers ---------------------------------------------


def get_intervals(interval_file: Intervals) -> list[tuple[str, int, int, str]]:
    """Read intervals from a BED file.

    ``track``/``browser``/comment/blank lines and lines with fewer than three
    columns are skipped.  A missing name column defaults to ``'.'``.

    Parameters
    ----------
    interval_file : str or Path
        Path to the BED file.

    Returns
    -------
    list of (str, int, int, str)
        ``(contig, start, stop, name)`` tuples.
    """
    intervals: list[tuple[str, int, int, str]] = []
    with open(interval_file, "r") as bed:
        for line in bed:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue

            contig = parts[0]
            start = int(parts[1])
            stop = int(parts[2])
            name = parts[3] if len(parts) > 3 else "."

            intervals.append((contig, start, stop, name))

    return intervals


def overlaps(
    contigs_1: NDArray,
    starts_1: NDArray,
    stops_1: NDArray,
    contigs_2: NDArray,
    starts_2: NDArray,
    stops_2: NDArray,
) -> NDArray:
    """Vectorized "does interval-set 1 overlap any interval in set 2?".

    Parameters
    ----------
    contigs_1, starts_1, stops_1 : ndarray
        Parallel arrays describing the query intervals.
    contigs_2, starts_2, stops_2 : ndarray
        Parallel arrays describing the reference intervals.

    Returns
    -------
    ndarray of bool
        Same shape as ``contigs_1``; ``True`` where that query interval
        overlaps at least one interval in set 2 (same contig required).
    """
    contigs_1 = contigs_1[:, np.newaxis]
    starts_1 = starts_1[:, np.newaxis]
    stops_1 = stops_1[:, np.newaxis]

    contigs_2 = contigs_2[np.newaxis]
    starts_2 = starts_2[np.newaxis]
    stops_2 = stops_2[np.newaxis]

    contig_blind_overlaps = np.logical_and(
        (starts_1 < stops_2), (stops_1 > starts_2)
    )
    in_same_contig = contigs_1 == contigs_2
    raw_overlaps = np.logical_and(contig_blind_overlaps, in_same_contig)
    return np.any(raw_overlaps, axis=1)


# -- k-mer helpers ----------------------------------------------------------


def gen_kmers(k: int, bases: str = "ACGT") -> list[str]:
    """Generate every length-``k`` k-mer over ``bases``.

    Parameters
    ----------
    k : int
        K-mer length (must be non-negative).
    bases : str, optional
        Alphabet, by default ``'ACGT'``.

    Returns
    -------
    list of str
        All ``len(bases) ** k`` k-mers, in lexicographic order.

    Raises
    ------
    ValueError
        If ``k`` is negative.
    """
    if k < 0:
        raise ValueError("k must be non-negative")
    return ["".join(p) for p in itertools.product(bases, repeat=k)]


@jit(nopython=True)
def _reverse_complement_numba(kmer_bytes: NDArray[np.uint8]) -> NDArray[np.uint8]:
    """Numba-accelerated reverse complement operating on ASCII bytes."""
    n = len(kmer_bytes)
    res = np.empty(n, dtype=np.uint8)
    for i in range(n):
        b = kmer_bytes[n - 1 - i]
        if b == 65 or b == 97:  # A / a
            res[i] = 84  # T
        elif b == 84 or b == 116:  # T / t
            res[i] = 65  # A
        elif b == 67 or b == 99:  # C / c
            res[i] = 71  # G
        elif b == 71 or b == 103:  # G / g
            res[i] = 67  # C
        else:
            res[i] = b  # keep (e.g. N)
    return res


def reverse_complement(kmer: str) -> str:
    """Return the reverse complement of a DNA string (``N`` preserved)."""
    kmer_bytes = np.frombuffer(kmer.encode("ascii"), dtype=np.uint8)
    rc_bytes = _reverse_complement_numba(kmer_bytes)
    return rc_bytes.tobytes().decode("ascii")
