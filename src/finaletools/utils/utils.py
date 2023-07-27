from __future__ import annotations
import time
import gzip
import tempfile as tf
from typing import Union, TextIO, Tuple, List, Generator
from sys import stderr, stdout

import numpy as np
from numpy.typing import NDArray
from numba import jit
import pysam
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
    in_region = np.logical_and(
        np.less(frag_array[:, 0], maximum),
        np.greater_equal(frag_array[:, 1], minimum)
        )
    filtered_frags = frag_array[in_region]
    return filtered_frags


def frag_generator(
    input_file: Union[str, pysam.AlignmentFile],
    contig: str,
    quality_threshold: int=30,
    start: int=None,
    stop: int=None,
    fraction_low: int=120,
    fraction_high: int=180,
    verbose: bool=False
) -> Generator[Tuple]:
    """
    Reads from BAM, SAM, or BED file and returns tuples containing
    contig (chromosome), start, and stop (end) for each fragment.

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
    frag_ends : Generator
        Generator that yields tuples containing the region covered by
        each fragment in input_file.
    """
    try:
        # check type of input and open if needed
        input_file_is_str = False   # file was opened in this context
        is_sam = False  # file is SAM/BAM, not tabix indexed
        if type(input_file) == str:   # path string
            input_file_is_str == True
            # check file type
            if (
                input_file.endswith('.sam')
                or input_file.endswith('.bam')
            ):
                is_sam = True
                sam_file = pysam.AlignmentFile(input_file, 'r')
            elif (
                input_file.endswith('frag.gz')
                or input_file.endswith('bed.gz')
                or input_file.endswith('frag.gz')
                or input_file.endswith('bed.gz')
            ):
                tbx = pysam.TabixFile(input_file, 'r')
        elif type(input_file) == pysam.AlignmentFile:
            is_sam = True
            sam_file = input_file
        elif type(input_file) == pysam.TabixFile:
            tbx = input_file
        else:
            raise TypeError(
                f'{type(input_file)} is invalid type for input_file.'
            )

        if is_sam:
            for read1 in sam_file.fetch(contig, start, stop):
                # Only select forward strand and filter out non-paired-end
                # reads and low-quality reads
                if (read1.is_read2
                    or low_quality_read_pairs(read1, quality_threshold)):
                    pass
                elif (
                    abs(read_length := read1.template_length) >= fraction_low
                    and abs(read_length) <= fraction_high
                ):
                    read_pos = read1.reference_start
                    read_mate_pos = read1.reference_start + read_length
                    if read_pos < read_mate_pos:
                        read_start, read_stop = read_pos, read_mate_pos
                    else:
                        read_start, read_stop = read_mate_pos, read_pos
                    read_on_plus = read1.is_forward
                    yield contig, read_start, read_stop, read_on_plus
        else:
            for line in tbx.fetch(
                contig, start, stop, parser=pysam.asTuple()
            ):
                read_start = int(line[1])
                read_stop = int(line[2])
                read_length = read_stop - read_start
                read_on_plus = '+' in line[5]
                if read_length >= fraction_low and read_length <= fraction_high:
                    yield contig, read_start, read_stop, read_on_plus
    finally:
        if input_file_is_str and is_sam:
            sam_file.close()
        elif input_file_is_str:
            tbx.close()


def frag_array(input_file: Union[str, pysam.AlignmentFile],
               contig: str,
               quality_threshold: int=30,
               start: int=None,
               stop: int=None,
               fraction_low: int=120,
               fraction_high: int=180,
               verbose: bool=False
               ) -> NDArray[np.int64]:
    """
    Reads from BAM, SAM, or BED file and returns a two column matrix
    with fragment start and stop positions.

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
    try:
        # check type of input and open if needed
        input_file_is_str = False   # file was opened in this context
        is_sam = False  # file is SAM/BAM, not tabix indexed
        if type(input_file) == str:   # path string
            input_file_is_str == True
            # check file type
            if (
                input_file.endswith('.sam')
                or input_file.endswith('.bam')
            ):
                is_sam = True
                sam_file = pysam.AlignmentFile(input_file, 'r')
            elif (
                input_file.endswith('frag.gz')
                or input_file.endswith('bed.gz')
                or input_file.endswith('frag.gz')
                or input_file.endswith('bed.gz')
            ):
                tbx = pysam.TabixFile(input_file, 'r')
        elif type(input_file) == pysam.AlignmentFile:
            is_sam = True
            sam_file = input_file
        elif type(input_file) == pysam.TabixFile:
            tbx = input_file
        else:
            raise TypeError(
                f'{type(input_file)} is invalid type for input_file.'
            )

        # based on file type, read into an array
        frag_ends = []
        if is_sam:
            for read1 in sam_file.fetch(contig, start, stop):
                # Only select forward strand and filter out non-paired-end
                # reads and low-quality reads
                if (read1.is_read2
                    or low_quality_read_pairs(read1, quality_threshold)):
                    pass
                # HACK: dealing with negative TLENs
                elif (
                    abs(read_length := read1.template_length) >= fraction_low
                    and abs(read_length) <= fraction_high
                ):
                    read_pos = read1.reference_start
                    read_mate_pos = read1.reference_start + read_length
                    read_start = min(read_pos, read_mate_pos)
                    read_stop = max(read_pos, read_mate_pos)
                    read_on_plus = read1.is_forward
                    frag_ends.append((read_start, read_stop, read_on_plus))
        else:
            for line in tbx.fetch(
                contig, start, stop, parser=pysam.asTuple()
            ):
                read_start = int(line[1])
                read_stop = int(line[2])
                read_length = read_stop - read_start
                read_on_plus = int('+' in line[5])
                if read_length >= fraction_low and read_length <= fraction_high:
                    frag_ends.append((read_start, read_stop, read_on_plus))
    finally:
        if input_file_is_str and is_sam:
            sam_file.close()
        elif input_file_is_str:
            tbx.close()

    # convert to ndarray
    frag_ends = np.array(frag_ends, dtype=np.int64)

    if frag_ends.ndim == 1:
        frag_ends = frag_ends.reshape((0, 2))

    assert frag_ends.ndim == 2, (f'frag_ends has dims {frag_ends.ndim} and '
                                 f'shape {frag_ends.shape}')
    assert (frag_ends.shape == (0, 2)
            or frag_ends.shape[1] == 2),('frag_ends has shape'
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


def genome2list(genome_file: str) -> list:
    """
    Reads a GENOME text file into a list of tuples (chrom, length)

    Parameters
    ----------
    genome_file : str
        String containing path to GENOME format file

    Returns
    _______
    chroms : str
        List of tuples containing chrom/contig names and lengths
    """
    chroms = []
    with open(genome_file) as file:
        for line in file:
            if line != '\n':
                chroms.append((
                    (contents:=line.split('\t'))[0],
                    int(contents[1])
                ))
    return chroms
