from __future__ import annotations
import time
import gzip
from typing import Union, Generator
from sys import stderr
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from numba import jit
import pysam


def chrom_sizes_to_list(
    chrom_sizes_file: Union[str, Path]) -> list[tuple[str][int]]:
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
    chrom_sizes_file: Union[str, Path]) -> list[tuple[str][int]]:
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


def _merge_overlapping_intervals(intervals):
    intervals.sort(key=lambda x: x[0])
    merged = []
    for interval in intervals:
        if not merged or interval[0] > merged[-1][1]:
            merged.append(interval)
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))
    return merged


def _reduce_overlaps_in_file(interval_file):
    intervals_dict = {}
    with open(interval_file, 'r') as file:
        for line in file:
            chrom, start, end = line.strip().split('\t')[:3]
            start, end = int(start), int(end)
            if chrom not in intervals_dict:
                intervals_dict[chrom] = []
            intervals_dict[chrom].append((start, end))

    reduced_intervals = {}
    for chrom, intervals in intervals_dict.items():
        reduced_intervals[chrom] = _merge_overlapping_intervals(intervals)
    return reduced_intervals    


def _convert_to_list(reduced_intervals):
    converted_intervals = {}
    for chrom, intervals in reduced_intervals.items():
        converted_intervals[chrom] = [[chrom, start, end] for start, end in intervals]
    return converted_intervals


def _merge_all_intervals(converted_intervals):
    all_intervals = []
    for intervals in converted_intervals.values():
        all_intervals.extend(intervals)
    return all_intervals


def frag_bam_to_bed(input_file: Union[str, pysam.AlignmentFile],
                    output_file: str,
                    contig: str=None,
                    quality_threshold: int=30,
                    verbose: bool=False):
    """
    Take paired-end reads from bam_file and write to a BED file.

    Parameters
    ----------
    input_file : pysam.AlignedFile or str
    output_file : str
    contig : str, optional
    quality_threshold : int, optional
    verbose : bool, optional
    """
    if (verbose):
        start_time = time.time()
        print('Opening file')

    sam_file = None
    try:
        # Open file or asign AlignmentFile to sam_file
        if (type(input_file) == pysam.AlignmentFile):
            sam_file = input_file
        elif (type(input_file) == str):
            sam_file = pysam.AlignmentFile(input_file)
        else:
            raise TypeError(
                ("bam_file should be an AlignmentFile or path string.")
                )

        # Open output file
        if output_file.endswith('.gz'):
            out = gzip.open(output_file, 'wt')
        else:
            out = open(output_file, 'w')

        # iterate through reads and send to BED
        for read1 in sam_file.fetch(contig=contig):
            # Only select forward strand and filter out non-paired-end
            # reads and low-quality reads
            if (_not_read1_or_low_quality(read1, quality_threshold)):
                pass
            else:
                out.write(
                    f'{read1.reference_name}\t{read1.reference_start}\t'
                    f'{read1.reference_start + read1.template_length}\n'
                    )
    except Exception as e:
        print("An error occurred:", str(e))

    finally:
        # Close everything when done
        if (type(input_file) == str):
            sam_file.close()
        out.close()

    if (verbose):
        end_time = time.time()
        print(f'frag_bam_to_bed took {end_time - start_time} s to complete',
              flush=True)


@jit(nopython=True)
def frags_in_region(frag_array: NDArray[np.int64],
                    minimum: int,
                    maximum: int) -> NDArray[np.int64]:
    """
    Takes an array of coordinates for ends of fragments and returns an
    array of fragments with coverage in the specified region. That is, a
    fragment is included if at least one base is in [minimum, maximum).

    Parameters
    ----------
    frag_array : ndarray
    minimum : int
    maximum : int

    Returns
    -------
    filtered_frags : ndarray
    """
    # Changed the code a bit to make it compatible with numba and not raise an error 
    starts = frag_array['start']
    stops = frag_array['stop']
    in_region = np.logical_and(np.less(starts,maximum), np.greater_equal(stops,minimum))
    filtered_frags = frag_array[in_region]
    return filtered_frags
def frag_generator(
    input_file: Union[str, pysam.AlignmentFile, pysam.TabixFile, Path],
    contig: str,
    quality_threshold: int=30,
    start: int=None,
    stop: int=None,
    fraction_low: int=120,
    fraction_high: int=180,
    intersect_policy: str="midpoint",
    verbose: bool=False
) -> Generator[tuple]:
    """
    Reads from BAM, SAM, or BED file and returns tuples containing
    contig (chromosome), start, stop (end), mapq, and strand for each fragment.
    Optionally may filter for mapq, size, and intersection with a region.

    Parameters
    ----------
    input_file : str or AlignmentFile
    contig : str
    quality_threshold : int, optional
    start : int, optional
    stop : int, optional
    fraction_low : int, optional
        Specifies lowest fragment length included in array. Default is
        120, equivalent to long fraction.
    fraction_high : int, optional
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
    frag_ends : Generator
        Generator that yields tuples containing the region covered by
        each fragment in input_file.
    """
    try:
        # check type of input and open if needed
        input_file_is_path = False
        if isinstance(input_file, str) or isinstance(input_file, Path):
            input_file_is_path == True
            # check file type
            if (
                str(input_file).endswith('.sam')
                or str(input_file).endswith('.bam')
                or str(input_file).endswith('.cram')
            ):
                is_sam = True
                sam_file = pysam.AlignmentFile(input_file, 'r')
            elif (
                str(input_file).endswith('frag.gz')
                or str(input_file).endswith('bed.gz')
            ):
                tbx = pysam.TabixFile(str(input_file), 'r')
                is_sam = False
            else:
                raise ValueError(
                    "Unaccepted interval file type. Only SAM, CRAM, BAM"
                    ", and Frag.gz files are accepted.")
        elif type(input_file) == pysam.AlignmentFile:
            input_file_is_path = False
            is_sam = True
            sam_file = input_file
        elif type(input_file) == pysam.TabixFile:
            input_file_is_path = False
            is_sam = False
            tbx = input_file
        else:
            raise TypeError(
                f'{type(input_file)} is invalid type for input_file.'
            )
            exit(1)
        
        # setting filter based on intersect policy
        if intersect_policy == 'midpoint':
            check_intersect = lambda r_start, r_stop, f_start, f_stop : (
                (r_start is None or ((f_start+f_stop)/2) >= r_start)
                and (r_stop is None or ((f_start+f_stop)/2) < r_stop)
            )
        elif intersect_policy == 'any':
            check_intersect = lambda r_start, r_stop, f_start, f_stop : (
                (r_start is None or f_stop > r_start)
                and (r_stop is None or f_start < r_stop)
            )
        else:
            raise ValueError(f'{intersect_policy} is not a valid policy')

        # Raise exception if start and stop specified but not contig
        if contig is None and not (start is None and stop is None):
            if contig is None and start==0 and stop is None:
                pass
            else:
                raise ValueError("contig should be specified if start or stop given.")

        if is_sam:
            for read in sam_file.fetch(contig, start, stop):
                # Only select read1 and filter out non-paired-end
                # reads and low-quality reads
                try:
                    if (low_quality_read_pairs(read, quality_threshold)
                        or read.is_read2):
                        pass
                    elif (
                        abs(frag_length := read.template_length) >= fraction_low
                        and abs(frag_length) <= fraction_high
                    ):
                        if read.template_length > 0:
                            f_start = read.reference_start
                            f_stop = read.reference_start + read.template_length
                            if (check_intersect(start, stop, f_start, f_stop)):
                                assert read.reference_start < read.reference_start + read.template_length, f"forward start {read.reference_start} after stop {read.reference_start + read.template_length} on chrom {read.reference_name} with read_is_forward {read.is_forward}."
                                yield (
                                    read.reference_name,
                                    read.reference_start,
                                    read.reference_start + read.template_length,
                                    read.mapping_quality,
                                    read.is_forward
                                )
                        elif read.template_length < 0:
                            f_start = read.reference_end + read.template_length
                            f_stop = read.reference_end
                            if (check_intersect(start, stop, f_start, f_stop)):
                                assert read.reference_end + read.template_length < read.reference_end, f"reverse start {read.reference_end + read.template_length} after stop {read.reference_end} on chrom {read.reference_name} with read_is_forward {read.is_forward}."
                                yield (
                                    read.reference_name,
                                    read.reference_end + read.template_length,
                                    read.reference_end,
                                    read.mapping_quality,
                                    read.is_forward 
                                )
                except TypeError as e:
                    stderr.writelines(["Type error encountered.\n",
                                       f"Fragment length: {frag_length}\n",
                                       f"fraction_low: {fraction_low}\n",
                                       f"fraction_high: {fraction_high}\n",
                                       "Skipping interval.\n",
                                       f"Error: {e}\n"])

        else:
            for line in tbx.fetch(
                contig, start, stop, parser=pysam.asTuple()
            ):
                read_start = int(line[1])
                read_stop = int(line[2])
                mapq = int(line[3])
                frag_length = read_stop - read_start
                read_on_plus = '+' in line[4]
                try:
                    if (frag_length >= fraction_low
                        and frag_length <= fraction_high
                        and mapq >= quality_threshold
                        ):
                        yield contig, read_start, read_stop, mapq, read_on_plus
                # HACK: for some reason read_length is sometimes None
                except TypeError:
                    continue

    finally:
        if input_file_is_path and is_sam:
            sam_file.close()
        elif input_file_is_path:
            tbx.close()


def frag_array(
    input_file: Union[str,pysam.AlignmentFile, pysam.TabixFile, Path],
    contig: str,
    quality_threshold: int=30,
    start: int=None,
    stop: int=None,
    fraction_low: int=120,
    fraction_high: int=180,
    intersect_policy: str="midpoint",
    verbose: bool=False
    ) -> NDArray:
    """
    Reads from BAM, SAM, or BED file and returns a three column matrix
    with fragment start and stop positions and strand.

    Parameters
    ----------
    input_file : str or AlignmentFile
    contig : str
    quality_threshold : int, optional
    start : int, optional
    stop : int, optional
    fraction_low : int, optional
        Specifies lowest fragment length included in array. Default is
        120, equivalent to long fraction.
    fraction_high : int, optional
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
            fraction_low,
            fraction_high,
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

    return (read.is_unmapped    # 0x4
            or read.is_secondary    # 0x100
            or (not read.is_paired) # not 0x1
            or read.mate_is_unmapped    # 0x8
            or read.is_duplicate    # 0x400
            or read.mapping_quality < min_mapq
            or read.is_qcfail   # 0x200
            or read.is_supplementary    # 0x800
            or (not read.is_proper_pair)   # not 0x2
            or (read.is_reverse and read.mate_is_reverse))   # -G 48


def _not_read1_or_low_quality(read: pysam.AlignedRead, min_mapq: int=30):
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


def _get_intervals(
    input_file: Union[str, pysam.AlignmentFile],
    interval_file: str,
    intersect_policy: str,
    quality_threshold: int,
    verbose: Union[bool, int]
) -> list[tuple[str, str, int, int, str, str, int]]:
    """
    Helper function to read intervals from bed file.
    Returns list of tuples:
    (input_file, chrom, start, stop, name, intersect_policy,
    quality_threshold, verbosity)
    """
    intervals = []  # list of inputs for single_coverage

    with open(interval_file) as bed:
        for line in bed:
            if not line.startswith('#'):
                if line != '':
                    contig, start, stop, *name = line.split()
                    start = int(start)
                    stop = int(stop)
                    name = name[0] if len(name) > 1 else '.'
                    interval = (
                        input_file,
                        contig,
                        start,
                        stop,
                        name,
                        intersect_policy,
                        quality_threshold,
                        verbose - 1 if verbose > 1 else 0
                    )
                    intervals.append(interval)
                else:
                    break
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