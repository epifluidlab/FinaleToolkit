from __future__ import annotations

import gzip
import os
import time
import warnings
from pathlib import Path
from sys import stderr
from typing import Generator

import numpy as np
import pysam
from numba import jit
from numpy.typing import NDArray

from finaletoolkit.utils.typing import ChromSizes

from ..io.alignment import AlignmentWrapper, Fragment
from ._comparison import _none_eq, _none_geq, _none_leq
from ._frag_generator import frag_generator


def chrom_sizes_to_list(
    chrom_sizes_file: ChromSizes
) -> list[tuple[str, int]]:
    """
    Reads chromosome names and sizes from a CHROMSIZE file into a list.

    Parameters
    ----------
    chrom_sizes_file: str or Path
        Tab-delimited file with column for chrom names and column for
        chrom sizes.
    
    Returns
    -------
    list of string, int tuples
        chrom names and sizes.
    """
    chrom_sizes = []
    with open(chrom_sizes_file, 'r') as file:
        for line in file:
            if line != '\n':
                chrom, size = line.strip().split('\t')
                chrom_sizes.append((chrom, int(size)))
    return chrom_sizes


def chrom_sizes_to_dict(
    chrom_sizes_file: ChromSizes) -> dict[str, int]:
    """
    Reads chromosome names and sizes from a CHROMSIZE file into a dict.

    Parameters
    ----------
    chrom_sizes_file: str or Path
        Tab-delimited file with column for chrom names and column for
        chrom sizes.
    
    Returns
    -------
    dict
        Chrom names are keys and values are int chrom sizes.
    """
    chrom_sizes = {}
    with open(chrom_sizes_file, 'r') as file:
        for line in file:
            if line != '\n':
                chrom, size = line.strip().split('\t')
                chrom_sizes[chrom] = int(size)
    return chrom_sizes


from ._intervals import (
    _merge_overlapping_intervals,
    _reduce_overlaps_in_file,
    _convert_to_list,
    _merge_all_intervals,
)


def frag_bam_to_bed(input_file: str | pysam.AlignmentFile,
                    output_file: str,
                    contig: str | None = None,
                    quality_threshold: int = 30,
                    verbose: bool = False,
                    reference_file: str | Path | None = None):
    """
    Take paired-end reads from bam_file and write to a BED file.

    Parameters
    ----------
    input_file : pysam.AlignedFile or str
    output_file : str
    contig : str, optional
    quality_threshold : int, optional
    verbose : bool, optional
    reference_file : str or Path, optional
        Reference genome file (required for CRAM).
    """
    if (verbose):
        start_time = time.time()
        print('Opening file')

    try:
        # Open output file
        if output_file.endswith('.gz'):
            out = gzip.open(output_file, 'wt')
        else:
            out = open(output_file, 'w')

        # Use AlignmentWrapper to handle file-type-specific logic
        with AlignmentWrapper(
            input_file, 
            reference_file=reference_file, 
            quality_threshold=quality_threshold
        ) as wrapper:
            # iterate through fragments and send to BED
            for frag in wrapper.fetch(contig=contig):
                out.write(f'{frag.contig}\t{frag.start}\t{frag.stop}\n')
    except Exception as e:
        logger.error("An error occurred during BAM to BED conversion: %s", str(e))
    finally:
        if 'out' in locals():
            out.close()

    if (verbose):
        end_time = time.time()
        print(f'frag_bam_to_bed took {end_time - start_time} s to complete',
              flush=True)


@jit(nopython=True)
def frags_in_region(frag_array: NDArray,
                    start: int,
                    stop: int) -> NDArray:
    """
    Takes an array of coordinates for ends of fragments and returns an
    array of fragments with coverage in the specified region. That is, a
    fragment is included if at least one base is in [minimum, maximum).

    Parameters
    ----------
    frag_array : ndarray
        Structured numpy array with 'start' and 'stop' fields.
    start : int
        Left-most coordinate of the region.
    stop : int
        Right-most coordinate of the region.

    Returns
    -------
    filtered_frags : ndarray
        The subset of fragments that overlap the region.
    """
    # Changed the code a bit to make it compatible with numba and not raise an error 
    starts = frag_array['start']
    stops = frag_array['stop']
    in_region = np.logical_and(np.less(starts, stop), np.greater_equal(stops, start))
    filtered_frags = frag_array[in_region]
    return filtered_frags


def frag_array(
    input_file: FragFile,
    contig: str,
    quality_threshold: int=30,
    start: int | None = None,
    stop: int | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    intersect_policy: str="midpoint",
    verbose: bool=False
    ) -> NDArray:
    """
    Reads from BAM, CRAM, or fragment file and returns a three column matrix
    with fragment start and stop positions and strand.

    Parameters
    ----------
    input_file : str or AlignmentFile
    contig : str
    quality_threshold : int, optional
    start : int, optional
    stop : int, optional
    min_length : int, optional
        Specifies lowest fragment length included in array. Default is
        120, equivalent to long fraction.
    max_length : int, optional
        Specifies highest fragment length included in array. Default is
        120, equivalent to long fraction.
        intersect_policy : str, optional
        Specifies what policy is used to include fragments in the
        given interval. Default is "midpoint". Policies include:
        - midpoint: the average of end coordinates of a fragment lies
        in the interval.
        - any: any part of the fragment is in the interval.
    verbose : bool, optional

    Returns
    -------
    frag_ends : NDArray
        'NDArray' with shape (n, 3) where column 1 contains fragment
        start position and column 2 contains fragment stop position, and
        column3 is 1 of on the + strand and is 0 if on the - strand.
        If no fragments exist in the specified minimum-maximum interval,
        the returned 'ndarray' will have a shape of (0, 3)
    """
    # use the frag_generator to create a list of intervals
    if verbose:
        stderr.write("frag_array: fetching fragments\n")

    frag_list = [
        (frag_start, frag_stop, strand)
        for _, frag_start, frag_stop, _, strand
        in frag_generator(
            input_file,
            contig,
            quality_threshold,
            start,
            stop,
            min_length,
            max_length,
            intersect_policy,
            verbose
        )
    ]
    if verbose:
        stderr.write("frag_array: converting to array\n")
    # convert to struct array
    frag_ends = np.array(frag_list, dtype=[("start", "i8"),("stop", "i8"),("strand", "?")])

    assert frag_ends.ndim == 1, (f'frag_ends has dims {frag_ends.ndim} and '
                                 f'shape {frag_ends.shape}')
    return frag_ends


def low_quality_read_pairs(read, min_mapq=30):
    """
    Return `True` if the sequenced read described in `read` is not a
    properly paired read with a Phred score exceeding `min_mapq`. Based
    on https://github.com/epifluidlab/cofragr/blob/master/python/frag_su
    mmary_in_intervals.py

    Equivalent to -F 3852 -f 3

    Parameters
    ----------
    read : pysam.AlignedSegment
        Sequenced read to check for quality, perfect pairing and if it
        is mapped.
    min_mapq : int, optional
        Minimum Phred score for map quality of read. Defaults to 30.

    Returns
    -------
    is_low_quality : bool
        True if read is low quality, unmapped, not properly paired.
    """

    if (
        read.is_unmapped    # 0x4
        or read.is_secondary    # 0x100
        or (not read.is_paired) # not 0x1
        or read.mate_is_unmapped    # 0x8
        or read.is_duplicate    # 0x400
        or read.mapping_quality < min_mapq
        or read.is_qcfail   # 0x200
        or read.is_supplementary    # 0x800
        or (not read.is_proper_pair)   # not 0x2
        or (read.is_reverse and read.mate_is_reverse)   # -G 48
    ):
        return True
    try:
        if read.has_tag("MQ") and read.get_tag("MQ") < min_mapq:
            return True
    except Exception:
        pass

    return False
    

def _not_read1_or_low_quality(read: pysam.AlignedSegment, min_mapq: int=30):
    """
    Return `True` if the sequenced read described in `read` is not read1
    and a properly paired read with a Phred score exceeding `min_mapq`.

    Parameters
    ----------
    read : pysam.AlignedSegment
        Sequenced read to check for quality, perfect pairing and if it
        is mapped.
    min_mapq : int, optional
        Minimum Phred score for map quality of read. Defaults to 30.

    Returns
    -------
    is_low_quality : bool
        True if read is either not read1, is low quality, is unmapped,
        or is not properly paired.
    """
    return (low_quality_read_pairs(read, min_mapq=min_mapq)
            or not read.is_read1)


def get_intervals(
    interval_file: Intervals
) -> list[tuple[str, int, int, str]]:
    """
    Helper function to read intervals from a BED file.

    Parameters
    ----------
    interval_file : str or Path
        Path to the BED file.

    Returns
    -------
    list of tuples
        Each tuple contains (contig, start, stop, name).
    """
    intervals = []
    with open(interval_file, 'r') as bed:
        for line in bed:
            if line.startswith(('#', 'track', 'browser')) or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
                
            contig = parts[0]
            start = int(parts[1])
            stop = int(parts[2])
            name = parts[3] if len(parts) > 3 else '.'
            
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
    """
    Function that performs vectorized computation of overlaps. Returns
    an array of same shape as contig_1 that is true if the intervals
    for set 1 each have any overlap with an interval in set 2.
    """
    contigs_1 = contigs_1[:, np.newaxis]
    starts_1 = starts_1[:, np.newaxis]
    stops_1 = stops_1[:, np.newaxis]

    contigs_2 = contigs_2[np.newaxis]
    starts_2 = starts_2[np.newaxis]
    stops_2 = stops_2[np.newaxis]

    contig_blind_overlaps = np.logical_and(
        (starts_1 < stops_2),
        (stops_1 > starts_2)
    )
    in_same_contig = contigs_1 == contigs_2
    raw_overlaps = np.logical_and(contig_blind_overlaps, in_same_contig)
    any_overlaps = np.any(raw_overlaps, axis=1)
    return any_overlaps

# None compatible comparison operators

import itertools
from typing import Generator

def gen_kmers(k: int, bases: str = 'ACGT') -> list[str]:
    """
    Generates all possible k-mers of length k using the given bases.
    Uses itertools.product for high-performance generation.

    Parameters
    ----------
    k : int
        Length of the k-mers to generate.
    bases : str, optional
        String containing the bases to use, by default 'ACGT'.

    Returns
    -------
    list[str]
        A list of all possible k-mer strings.
    """
    if k < 0:
        raise ValueError("k must be non-negative")
    
    return [''.join(p) for p in itertools.product(bases, repeat=k)]

@jit(nopython=True)
def _reverse_complement_numba(kmer_bytes: NDArray[np.uint8]) -> NDArray[np.uint8]:
    """Numba-optimized reverse complement using bytes."""
    n = len(kmer_bytes)
    res = np.empty(n, dtype=np.uint8)
    for i in range(n):
        b = kmer_bytes[n - 1 - i]
        if b == 65 or b == 97: # A or a
            res[i] = 84 # T
        elif b == 84 or b == 116: # T or t
            res[i] = 65 # A
        elif b == 67 or b == 99: # C or c
            res[i] = 71 # G
        elif b == 71 or b == 103: # G or g
            res[i] = 67 # C
        else:
            res[i] = b # Keep as is (e.g. N)
    return res

def reverse_complement(kmer: str) -> str:
    """Reverse complement a k-mer string using Numba optimization."""
    kmer_bytes = np.frombuffer(kmer.encode('ascii'), dtype=np.uint8)
    rc_bytes = _reverse_complement_numba(kmer_bytes)
    return rc_bytes.tobytes().decode('ascii')
