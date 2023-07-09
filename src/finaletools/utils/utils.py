from __future__ import annotations
import time
import gzip
import tempfile as tf
from typing import Union, TextIO, Tuple, List
from sys import stderr, stdout

import numpy as np
from numba import jit
import pysam
import pybedtools
from tqdm import tqdm


def _get_contigs(
        input_file: Union[str, pysam.AlignmentFile],
        verbose: bool=False
    ) -> list:
    """
    Retrieves contigs from input_file and returns lists of contig names
    and lengths
    """

    input_is_file = False
    try:
        # handling input types
        if (type(input_file) == pysam.AlignmentFile):
            sam_file = input_file
        elif input_file.endswith('bam'):
            input_is_file = True
            if (verbose):
                stderr.write(f'Opening {input_file}\n')
            sam_file = pysam.AlignmentFile(input_file)
        else:
            raise ValueError(
                'Invalid input_file type. Only BAM or SAM files are allowed.'
            )
        contigs = sam_file.references
        lengths = sam_file.lengths
    finally:
        if input_is_file:
            sam_file.close()

    return zip(contigs, lengths)


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
def frags_in_region(frag_array: np.ndarray,
                   minimum: int,
                   maximum: int) -> np.ndarray:
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
    in_region = np.logical_and(
        np.less(frag_array[:, 0], maximum),
        np.greater_equal(frag_array[:, 1], minimum)
        )
    filtered_frags = frag_array[in_region]
    return filtered_frags


def _sam_frag_array(sam_file: pysam.AlignmentFile,
                    contig: str,
                    has_min_max: bool,
                    quality_threshold:int=30,
                    minimum: int=None,
                    maximum: int=None,
                    fraction_low: int=120,
                    fraction_high: int=180,
                    verbose: bool=False):
    frag_ends = []
    count = sam_file.count(
        contig=contig,
        start=minimum,
        stop=maximum) if verbose else None
    if (has_min_max):
        for read1 in (tqdm(
            sam_file.fetch(contig=contig, start=minimum, stop=maximum),
            total=count)
            if verbose
            else sam_file.fetch(contig=contig)
        ):
            # Only select forward strand and filter out non-paired-end
            # reads and low-quality reads
            if (read1.is_read2
                or low_quality_read_pairs(read1, quality_threshold)):
                pass
            else:
                read_length = read1.template_length
                read_start = read1.reference_start
                read_stop = read1.reference_start + read_length
                if ((read_length >= fraction_low)
                    and (read_length <= fraction_high)):
                    frag_ends.append((read_start, read_stop))
    elif (contig is None):
        for read1 in (tqdm(sam_file, total=count)
                      if verbose
                      else sam_file):
            # Only select forward strand and filter out non-paired-end
            # reads and low-quality reads
            if (read1.is_read2
                or low_quality_read_pairs(read1, quality_threshold)):
                pass
            else:
                frag_ends.append(
                    (read1.reference_start,
                     read1.reference_start + read1.template_length)
                     )
    else:
        for read1 in (tqdm(sam_file.fetch(contig=contig), total=count)
                      if verbose
                      else sam_file.fetch(contig=contig)):
            # Only select forward strand and filter out non-paired-end
            # reads and low-quality reads
            if (read1.is_read2
                or low_quality_read_pairs(read1, quality_threshold)):
                pass
            else:
                frag_ends.append(
                    (read1.reference_start,
                     read1.reference_start + read1.template_length)
                     )
    return frag_ends


@jit(forceobj=True)
def _bed_frag_array(bed_file: TextIO,
                    contig: str,
                    has_min_max: bool,
                    quality_threshold:int=30,
                    minimum: int=None,
                    maximum: int=None,
                    fraction_low: int=120,
                    fraction_high: int=180,
                    verbose: bool=False):
    frag_ends = []
    if (has_min_max):
        for line in bed_file:
            frag_info = line.split('\t')
            read_start = int(frag_info[1])
            read_stop = int(frag_info[2])
            read_length = read_stop - read_start
            if ((frag_info[0] == contig)
                and ((read_stop >= minimum)
                     and (read_start < maximum)
                     and (read_length >= fraction_low)
                     and (read_length <= fraction_high)
                     )
                ):
                frag_ends.append((read_start, read_stop))
    else:
        for line in bed_file:
            frag_info = line.split('\t')
            if (frag_info[0] == contig):
                frag_ends.append((int(frag_info[1]), int(frag_info[2])))
    return frag_ends


def frag_array(input_file: Union[str, pysam.AlignmentFile],
               contig: str,
               quality_threshold: int=30,
               minimum: int=None,
               maximum: int=None,
               fraction_low: int=120,
               fraction_high: int=180,
               verbose: bool=False
               ) -> np.ndarray:
    """
    Reads from BAM, SAM, or BED file and returns a two column matrix
    with fragment start and stop positions.

    Parameters
    ----------
    input_file : str or AlignmentFile
    contig : str
    quality_threshold : int, optional
    minimum : int, optional
    maximum : int, optional
    fraction_low : int, optional
        Specifies lowest fragment length included in array. Default is
        120, equivalent to long fraction.
    fraction_high : int, optional
        Specifies highest fragment length included in array. Default is
        120, equivalent to long fraction.
    verbose : bool, optional

    Returns
    -------
    frag_ends : ndarray
        'ndarray' with shape (n, 2) where column 1 contains fragment
        start positions and column 2 contains fragment stop positions.
        If no fragments exist in the specified minimum-maximum interval,
        the returned 'ndarray' will have a shape of (0, 2)
    """
    # boolean flag indicating whether or not a minimum and maximum
    # location for one fragment end is specified.
    has_min_max = (minimum is not None) and (maximum is not None)
    if (minimum is None) != (maximum is None):
        raise ValueError(
            'Both minimum and maximum must either be present or absent.'
            )

    # input_file is AlignmentFile
    if (type(input_file) == pysam.AlignmentFile):
        sam_file = input_file
        frag_ends = _sam_frag_array(sam_file,
                                    contig,
                                    has_min_max,
                                    quality_threshold=quality_threshold,
                                    minimum=minimum, maximum=maximum,
                                    fraction_low=fraction_low,
                                    fraction_high=fraction_high,
                                    verbose=verbose)

    # input_file is a path string
    elif (type(input_file) == str):
        # SAM file
        if (input_file.endswith('.bam')
            or input_file.endswith('.sam')
            or input_file == '-'
        ):
            with pysam.AlignmentFile(input_file, 'r') as sam_file:
                frag_ends = _sam_frag_array(
                sam_file,
                contig,
                has_min_max,
                quality_threshold=quality_threshold,
                minimum=minimum,
                maximum=maximum,
                fraction_low=fraction_low,
                fraction_high=fraction_high,
                verbose=verbose)
        # BAM file
        elif (input_file.endswith('.bam')):
            with pysam.AlignmentFile(input_file, 'rb') as sam_file:
                frag_ends = _sam_frag_array(
                sam_file,
                contig,
                has_min_max,
                quality_threshold=quality_threshold,
                minimum=minimum,
                maximum=maximum,
                fraction_low=fraction_low,
                fraction_high=fraction_high,
                verbose=verbose)
        # BAM from stdin
        elif (input_file == '-'):
            with pysam.AlignmentFile(input_file, 'rb') as sam_file:
                frag_ends = _sam_frag_array(
                sam_file,
                contig=None
                False,
                quality_threshold=quality_threshold,
                fraction_low=fraction_low,
                fraction_high=fraction_high,
                verbose=verbose)
        # BED file
        elif (input_file.endswith('.bed')):
            with open(input_file, 'rt') as bed_file:
                frag_ends = _bed_frag_array(
                bed_file,
                contig,
                has_min_max,
                quality_threshold=quality_threshold,
                minimum=minimum,
                maximum=maximum,
                fraction_low=fraction_low,
                fraction_high=fraction_high,
                verbose=verbose)
        # BED.gz file
        elif (input_file.endswith('.bed.gz')):
            with gzip.open(input_file, 'rt') as bed_file:
                frag_ends = _bed_frag_array(
                    bed_file,
                    contig,
                    has_min_max,
                    quality_threshold=quality_threshold,
                    minimum=minimum,
                    maximum=maximum,
                    fraction_low=fraction_low,
                    fraction_high=fraction_high,
                    verbose=verbose)
        else:
            raise ValueError(
                'input_file can only have suffixes .bam, .sam, .bed, or '
                '.bed.gz'
                )

    else:
        raise TypeError(
            f'input_file is unsupported type "{type(input_file)}". Input_file '
            'should be a pysam.AlignmentFile or a string containing the path '
            'to a SAM or BAM file.'
            )

    # convert to ndarray
    frag_ends = np.array(frag_ends)

    if frag_ends.ndim == 1:
        frag_ends = frag_ends.reshape((0, 2))

    assert frag_ends.ndim == 2, (f'frag_ends has dims {frag_ends.ndim} and '
                                 f'shape {frag_ends.shape}')
    assert (frag_ends.shape == (0, 2)
            or frag_ends.shape[1] == 2), ('frag_ends has shape'
                                          f'{frag_ends.shape}')
    return frag_ends


def low_quality_read_pairs(read, min_mapq=30):
    """
    Return `True` if the sequenced read described in `read` is not a
    properly paired read with a Phred score exceeding `min_mapq`. Based
    on https://github.com/epifluidlab/cofragr/blob/master/python/frag_su
    mmary_in_intervals.py

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

    return (read.is_unmapped
            or read.is_secondary
            or (not read.is_paired)
            or read.mate_is_unmapped
            or read.is_duplicate
            or read.mapping_quality < min_mapq
            or read.is_qcfail
            or read.is_supplementary
            or (not read.is_proper_pair)
            or read.reference_name != read.next_reference_name)


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
    quality_threshold: int,
    verbose: Union[bool, int]
) -> list:
    """Helper function to read intervals from bed file."""
    intervals = []  # list of inputs for single_coverage

    with open(interval_file) as bed:
        for line in bed:
            if ~line.startswith('#'):
                if line != '':
                    contents = line.split()
                    contig = contents[0].strip()
                    start = int(contents[1])
                    stop = int(contents[2])
                    name = contents[3] if len(contents) > 3 else '.'
                    interval = (
                        input_file,
                        contig,
                        start,
                        stop,
                        name,
                        quality_threshold,
                        verbose - 1 if verbose > 1 else 0
                    )
                    intervals.append(interval)
                else:
                    break
    return intervals
