from __future__ import annotations
import gzip
import time
import pysam
import numpy as np
from numba import jit
from tqdm import tqdm
import tempfile as tf
from sys import stderr, stdout
from numpy.typing import NDArray
from typing import Union, TextIO, Tuple, List, Generator

def _parse_chrom_sizes(chrom_sizes_file):
    chrom_sizes = []
    with open(chrom_sizes_file, 'r') as file:
        for line in file:
            chrom, size = line.strip().split('\t')
            chrom_sizes.append((chrom, int(size)))
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

def get_contig_lengths(
    input_file: Union[str, pysam.AlignmentFile],
    verbose: bool
) -> List[Tuple[str, int]]:
    
    input_is_file = False
    sam_file = None
    
    try:
        if isinstance(input_file, pysam.AlignmentFile):
            sam_file = input_file
        elif isinstance(input_file, str) and input_file.endswith(('bam', 'sam')):
            input_is_file = True
            if verbose:
                stderr.write(f'Opening {input_file}\n')
            sam_file = pysam.AlignmentFile(input_file)
        else:
            raise ValueError('Invalid input_file type. Only BAM or SAM files are allowed.')

        contigs = sam_file.references
        lengths = sam_file.lengths

        contig_lengths = list(zip(contigs, lengths))

    finally:
        if input_is_file and sam_file is not None:
            sam_file.close()

    return contig_lengths

def chromsizes2list(chrom_sizes_file: str) -> list:
    chroms = []
    with open(chrom_sizes_file) as file:
        for line in file:
            if line != '\n':
                chroms.append((
                    (contents:=line.split('\t'))[0],
                    int(contents[1])
                ))
    return chroms

@jit(nopython=True)
def frags_in_region(frag_array: NDArray[np.int64],
                    minimum: int,
                    maximum: int) -> NDArray[np.int64]:
    in_region = np.logical_and(
        np.less(frag_array[:, 0], maximum),
        np.greater_equal(frag_array[:, 1], minimum)
    )
    filtered_frags = frag_array[in_region]
    return filtered_frags

def _get_file_handle(input_file: Union[str, pysam.AlignmentFile, pysam.TabixFile], verbose: bool):
    file_handle = None
    file_type = None

    if isinstance(input_file, pysam.AlignmentFile):
        file_handle = input_file
        file_type = 'sam'
    elif isinstance(input_file, pysam.TabixFile):
        file_handle = input_file
        file_type = 'tabix'
    elif isinstance(input_file, str):
        if input_file.endswith(('.sam', '.bam', '.cram')):
            file_handle = pysam.AlignmentFile(input_file, 'r')
            file_type = 'sam'
        elif input_file.endswith(('.frag.gz', '.bed.gz')):
            file_handle = pysam.TabixFile(input_file, 'r')
            file_type = 'tabix'
        else:
            raise ValueError('Invalid file format. Supported formats are SAM/BAM/CRAM or indexed bed/frag.gz files.')

        if verbose:
            stderr.write(f'Opening {input_file}\n')

    return file_handle, file_type

def _close_file_handle(file_handle, file_type):
    if file_type == 'sam' or file_type == 'tabix':
        file_handle.close()

def _check_intersect_policy(intersect_policy):
    if intersect_policy == 'midpoint':
        return lambda r_start, r_stop, f_start, f_stop: ((r_start is None or ((f_start+f_stop)/2) >= r_start)
                                                        and (r_stop is None or ((f_start+f_stop)/2) < r_stop))
    elif intersect_policy == 'any':
        return lambda r_start, r_stop, f_start, f_stop: ((r_start is None or f_stop > r_start)
                                                        and (r_stop is None or f_start < r_stop))
    else:
        raise ValueError(f'{intersect_policy} is not a valid policy. Please use "midpoint" or "any".')

def frag_generator(input_file: Union[str, pysam.AlignmentFile, pysam.TabixFile],
                   contig: str,
                   quality_threshold: int,
                   start: int,
                   stop: int,
                   min: int,
                   max: int,
                   intersect_policy: str,
                   verbose: bool):

    sam_file, file_type = _get_file_handle(input_file, verbose)
    try:
        if file_type == 'sam':
            check_intersect = _check_intersect_policy(intersect_policy)
            for read in sam_file.fetch(contig, start, stop):
                if (low_quality_read_pairs(read, quality_threshold) or read.is_read2):
                    continue
                frag_length = abs(read.template_length)
                if (min is None or min <= frag_length) and (max is None or frag_length <= max):
                    f_start, f_stop = (read.reference_start, read.reference_start + read.template_length) if read.template_length > 0 else (
                        read.reference_end + read.template_length, read.reference_end)
                    if check_intersect(start, stop, f_start, f_stop):
                        if f_start !=f_stop:
                            yield (read.reference_name, f_start, f_stop, read.mapping_quality, read.is_forward)
        elif file_type == 'tabix':
            check_intersect = _check_intersect_policy(intersect_policy)
            for line in sam_file.fetch(contig, start, stop, parser=pysam.asTuple()):
                read_start, read_stop, mapq = int(line[1]), int(line[2]), int(line[3])
                frag_length = read_stop - read_start
                if (check_intersect(start, stop, read_start, read_stop)):
                    if (min is None or min <= frag_length) and (max is None or frag_length <= max) and mapq >= quality_threshold:
                        yield (contig, read_start, read_stop, mapq, '+' in line[4])
    finally:
        _close_file_handle(sam_file, file_type)

def frag_array(input_file: Union[str, pysam.AlignmentFile],
               contig: str,
               quality_threshold: int,
               start: int,
               stop: int,
               min: int,
               max: int,
               intersect_policy: str,
               verbose: bool) -> np.ndarray:

    if verbose:
        stderr.write("frag_array: Fetching fragments.\n")

    frag_list = [(frag_start, frag_stop, strand) for _, frag_start, frag_stop, _, strand in frag_generator(input_file, contig, quality_threshold, start, stop, min, max, intersect_policy, verbose)]

    if verbose:
        stderr.write("frag_array: Converting to array.\n")

    frag_ends = np.array(frag_list, dtype=[("start", "i8"), ("stop", "i8"), ("strand", "?")])

    assert frag_ends.ndim == 1, f'frag_array has dims {frag_ends.ndim} and shape {frag_ends.shape}'
    return frag_ends

def low_quality_read_pairs(read, min_mapq):
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

def _not_read1_or_low_quality(read: pysam.AlignedRead, min_mapq: int):
    return (low_quality_read_pairs(read, min_mapq=min_mapq)
            or not read.is_read1)

def _get_intervals(
    interval_file: str,
) -> List[Tuple[str, str, int, int, str, str, int]]:
    
    intervals = []  # List of inputs for single_coverage

    with open(interval_file) as bed:
        for line in bed:
            if not line.startswith('#'):
                if line.strip():  # Check if line is not empty
                    contig, start, stop, *name = line.split()
                    start = int(start)
                    stop = int(stop)
                    name = name[0] if name else '.'
                    interval = (
                        contig,
                        start,
                        stop,
                        name
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
