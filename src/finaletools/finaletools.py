"""
Author: James Li
Created: 6/2/23
PI: Yaping Liu

Description:
Python script to calculate fragment features given a BAM file.

"""
# TODO: typing annotations for all functions

from __future__ import annotations
import pysam
import argparse
import gzip
import numpy as np
import time
import tempfile as tf
from numba import jit
from tqdm import tqdm
from multiprocessing.pool import Pool
from typing import Union, TextIO, BinaryIO


def frag_bam_to_bed(input_file,
                    output_file,
                    contig=None,
                    quality_threshold=15,
                    verbose=False):
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
            if (not_read1_or_low_quality(read1, quality_threshold)):
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
        print(f'frag_bam_to_bed took {end_time - start_time} s to complete')


def frags_in_region(frag_array: np.ndarray[int, int],
                   minimum: int,
                   maximum: int) -> np.ndarray[int, int]:
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
                    quality_threshold:int=15,
                    minimum: int=None,
                    maximum: int=None,
                    fraction_low: int=120,
                    fraction_high: int=180,
                    verbose: bool=False):
    frag_ends = []
    count = sam_file.count(contig=contig) if verbose else None
    if (has_min_max):
        for read1 in (tqdm(sam_file.fetch(contig=contig), total=count)
                      if verbose
                      else sam_file.fetch(contig=contig)):
            # Only select forward strand and filter out non-paired-end
            # reads and low-quality reads
            if (read1.is_read2
                or low_quality_read_pairs(read1, quality_threshold)):
                pass
            else:
                read_length = read1.template_length
                read_start = read1.reference_start
                read_stop = read1.reference_start + read_length
                if ((read_stop >= minimum)
                    and (read_start < maximum)
                    and (read_length >= fraction_low)
                    and (read_length <= fraction_high)):
                    frag_ends.append((read_start, read_stop))
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


@jit
def _bed_frag_array(bed_file: TextIO,
                    contig: str,
                    has_min_max: bool,
                    quality_threshold:int=15,
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
               quality_threshold: int=15,
               minimum: int=None,
               maximum: int=None,
               fraction_low: int=120,
               fraction_high: int=180,
               verbose: bool=False
               ) -> np.ndarray[int, int]:
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

    # input_file is AllignmentFile
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
        # BAM or SAM file
        if (input_file.endswith('.bam') or input_file.endswith('.sam')):
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


def low_quality_read_pairs(read, min_mapq=15):
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
        Minimum Phred score for map quality of read. Defaults to 15.

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


def not_read1_or_low_quality(read: pysam.AlignedRead, min_mapq: int=15):
    """
    Return `True` if the sequenced read described in `read` is not read1
    and a properly paired read with a Phred score exceeding `min_mapq`.

    Parameters
    ----------
    read : pysam.AlignedSegment
        Sequenced read to check for quality, perfect pairing and if it
        is mapped.
    min_mapq : int, optional
        Minimum Phred score for map quality of read. Defaults to 15.

    Returns
    -------
    is_low_quality : bool
        True if read is either not read1, is low quality, is unmapped,
        or is not properly paired.
    """
    return (low_quality_read_pairs(read, min_mapq=min_mapq)
            or not read.is_read1)


def frag_length(input_file: Union[str, pysam.AlignedSegment],
                contig: str=None,
                output_file: str=None, workers: int=1,
                quality_threshold: int=15,
                verbose: bool=False
                ) -> np.ndarray[int]:
    """
    Return `np.ndarray` containing lengths of fragments in `input_file`
    that are above the quality threshold and are proper-paired reads.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, SAM, or CRAM file containing paired-end fragment reads or
        its path. `AlignmentFile` must be opened in read mode.
    contig : string, optional
    output_file : string, optional
    quality_threshold : int, optional
    verbose : bool, optional

    Returns
    -------
    lengths : numpy.ndarray
        `ndarray` of fragment lengths from file and contig if
        specified.
    """
    if (verbose):
        start_time = time.time()
        print("Finding frag lengths.")

    lengths = []    # list of fragment lengths
    if (type(input_file) == pysam.AlignmentFile):
        sam_file = input_file
        if (verbose):
            print('Counting reads')
        count = sam_file.count(contig=contig) if verbose else None
        if (verbose):
            print(f'{count} reads counted')
        # Iterating on each read in file in specified contig/chromosome
        for read1 in (tqdm(sam_file.fetch(contig=contig), total=count)
                      if verbose
                      else sam_file.fetch(contig=contig)):
            # Only select forward strand and filter out non-paired-end
            # reads and low-quality reads
            if (not_read1_or_low_quality(read1, quality_threshold)):
                pass
            else:
                # append length of fragment to list
                lengths.append(abs(read1.template_length))
    else:
        if (verbose):
            print(f'Opening {input_file}')
        with pysam.AlignmentFile(input_file) as sam_file:   # Import
            if (verbose):
                print('Counting reads')
            count = sam_file.count(contig=contig) if verbose else None
            if (verbose):
                print(f'{count} reads counted')
            # Iterating on each read in file in specified
            # contig/chromosome
            for read1 in (tqdm(sam_file.fetch(contig=contig), total=count)
                          if verbose
                          else sam_file.fetch(contig=contig)):
                # Only select forward strand and filter out
                # non-paired-end reads and low-quality reads
                if (not_read1_or_low_quality(read1, quality_threshold)):
                    pass
                else:
                    # append length of fragment to list
                    lengths.append(abs(read1.template_length))

    # convert to array
    lengths = np.array(lengths)

    # check if output specified
    if (type(output_file) == str):
        if output_file.endswith(".bin"): # binary file
            with open(output_file, 'wt') as out:
                lengths.tofile(out)
        else:   # unaccepted file type
            raise ValueError(
                'output_file can only have suffixes .wig or .wig.gz.'
                )

    elif (output_file is not None):
        raise TypeError(
            f'output_file is unsupported type "{type(input_file)}". '
            'output_file should be a string specifying the path of the file '
            'to write output scores to.'
            )

    if (verbose):
        end_time = time.time()
        print(f'frag_length took {end_time - start_time} s to complete')

    return lengths


def frag_center_coverage(input_file,
                         contig,
                         start,
                         stop,
                         output_file=None,
                         quality_threshold=15,
                         verbose=False):
    """
    Return estimated fragment coverage over specified `contig` and
    region of`input_file`. Uses an algorithm where the midpoints of
    fragments are calculated and coverage is tabulated from the
    midpoints that fall into the specified region. Not suitable for
    fragments of size approaching region size.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, SAM, or CRAM file containing paired-end fragment reads or
        its path. `AlignmentFile` must be opened in read mode.
    contig : string
    start : int
    stop : int
    output_file : string, optional
    quality_threshold : int, optional
    verbose : bool, optional

    Returns
    -------
    coverage : int
        Fragment coverage over contig and region.
    """
    # TODO: determine if reference (like as found in pysam) is necessary
    # TODO: consider including region string (like in pysam)

    if (verbose):
        start_time = time.time()

    # initializing variable for coverage tuple outside of with statement
    coverage = 0

    if (type(input_file) == pysam.AlignmentFile):
        sam_file = input_file
        # Iterating on each read in file in specified contig/chromosome
        for read1 in sam_file.fetch(contig=contig):
            # Only select forward strand and filter out non-paired-end
            # reads and low-quality reads
            if (not_read1_or_low_quality(read1, quality_threshold)):
                pass
            else:
                # calculate mid-point of fragment
                center = read1.reference_start + read1.template_length // 2
                if ((center >= start) and (center < stop)):
                    coverage += 1
    else:
        with pysam.AlignmentFile(input_file, 'r') as sam_file:
            # Iterating on each read in file in
            # specified contig/chromosome
            for read1 in sam_file.fetch(contig=contig):
                # Only select forward strand and filter out
                # non-paired-end reads and low-quality reads
                if not_read1_or_low_quality(read1, quality_threshold):
                    pass
                else:
                    # calculate mid-point of fragment
                    center = read1.reference_start + read1.template_length // 2
                    if ((center >= start) and (center < stop)):
                        coverage += 1
    if (verbose):
        end_time = time.time()
        print(f'frag_coverage took {end_time - start_time} s to complete')

    return coverage


@jit(nopython=True)
def _single_wps(window_start: int,
                window_stop: int,
                window_position: int,
                frag_ends: np.ndarray[int, int]
                ) -> tuple[int, int]:
    # count number of totally spanning fragments
    is_spanning = ((frag_ends[:, 0] < window_start)
                   * (frag_ends[:, 1] > window_stop))
    num_spanning = np.sum(is_spanning)

    # count number of fragments with end in window
    is_start_in = ((frag_ends[:, 0] >= window_start)
                   * (frag_ends[:, 0] <= window_stop))
    is_stop_in = ((frag_ends[:, 1] >= window_start)
                  * (frag_ends[:, 1] <= window_stop))
    is_end_in = np.logical_or(is_start_in, is_stop_in)
    num_end_in = np.sum(is_end_in)

    # calculate wps and return
    return (window_position, num_spanning - num_end_in)


@jit(nopython=True)
def _vectorized_wps(frag_ends, window_starts, window_stops):
    """
    Unused helper function for vectorization
    """


    w_starts = np.column_stack(window_starts)
    w_stops = np.column_stack(window_stops)
    frag_starts = np.row_stack(frag_ends[:, 0])
    frag_stops = np.row_stack(frag_ends[:, 1])

    is_spanning = np.logical_and(
            np.less_equal(frag_starts, w_starts),
            np.greater_equal(frag_stops, w_stops))

    n_spanning = np.sum(is_spanning, axis=0)

    start_in = np.logical_and(
        np.less(frag_starts, w_starts),
        np.greater_equal(frag_starts, w_stops))

    stop_in = np.logical_and(
        np.less(frag_stops, w_starts),
        np.greater_equal(frag_stops, w_stops))

    end_in = np.logical_or(start_in, stop_in)

    n_end_in = np.sum(end_in, axis=0)

    scores = n_spanning - n_end_in

    return scores


@jit(nopython=True)
def _wps_loop(frag_ends: np.ndarray[int],
              start: int,
              stop: int,
              window_size: int):
    # array to store positions and scores
    scores = np.zeros((stop-start, 2))
    window_centers = np.arange(start, stop, dtype=np.int64)
    scores[:, 0] = window_centers
    window_starts = np.zeros(stop-start)
    window_stops = np.zeros(stop-start)
    np.rint(window_centers - window_size * 0.5, window_starts)
    np.rint(window_centers + window_size * 0.5 - 1, window_stops)
    # inclusive

    for i in range(stop-start):
        scores[i, :] = _single_wps(
            window_starts[i],
            window_stops[i],
            window_centers[i],
            frag_ends)

    return scores


def wps(input_file: Union[str, pysam.AlignmentFile],
        contig: str,
        start: Union[int, str],
        stop: Union[int, str],
        output_file: str=None,
        window_size: int=120,
        fraction_low: int=120,
        fraction_high: int=180,
        quality_threshold: int=15,
        verbose: Union[bool, int]=0
        ) -> np.ndarray[int, int]:
    """
    Return Windowed Protection Scores as specified in Snyder et al
    (2016) over a region [start,stop).

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM or SAM file containing paired-end fragment reads or its
        path. `AlignmentFile` must be opened in read mode.
    contig : str
    start : int
    stop : int
    output_file : string, optional
    window_size : int, optional
        Size of window to calculate WPS. Default is k = 120, equivalent
        to L-WPS.
    fraction_low : int, optional
        Specifies lowest fragment length included in calculation.
        Default is 120, equivalent to long fraction.
    fraction_high : int, optional
        Specifies highest fragment length included in calculation.
        Default is 120, equivalent to long fraction.
    quality_threshold : int, optional
    workers : int, optional
    verbose : bool, optional

    Returns
    -------
    scores : numpy.ndarray
        np array of shape (n, 2) where column 1 is the coordinate and
        column 2 is the score and n is the number of coordinates in
        region [start,stop)
    """

    if (verbose):
        start_time = time.time()
        print("Reading fragments")

    # set start and stop to ints
    start = int(start)
    stop = int(stop)

    # set minimum and maximum values for fragments. These extend farther
    # than needed
    minimum = round(start - window_size)
    maximum = round(stop + window_size)

    # read fragments from file
    frag_ends = frag_array(input_file,
                           contig,
                           quality_threshold,
                           minimum=minimum,
                           maximum=maximum,
                           fraction_low=fraction_low,
                           fraction_high=fraction_high,
                           verbose=(verbose>=2))

    if (verbose):
        print("Done reading fragments, preparing for WPS calculation.")
    # check if no fragments exist on this interval
    if (frag_ends.shape == (0, 2)):

        scores = np.zeros((stop-start, 2))
        scores[:, 0] = np.arange(start, stop, dtype=int)
    else:
        scores = _wps_loop(frag_ends, start, stop, window_size)


    # TODO: consider switch-case statements and determine if they
    # shouldn't be used for backwards compatability
    if (type(output_file) == str):   # check if output specified

        if (verbose):
            print('Writing to output file.')

        if output_file.endswith(".wig.gz"): # zipped wiggle
            with gzip.open(output_file, 'wt') as out:
                # declaration line
                out.write(
                    f'fixedStep\tchrom={contig}\tstart={start}\t'
                    f'step={1}\tspan={stop-start}\n'
                    )
                for score in scores[:, 1]:
                    out.write(f'{score}\n')

        elif output_file.endswith(".wig"):  # wiggle
            with open(output_file, 'wt') as out:
                # declaration line
                out.write(
                    f'fixedStep\tchrom={contig}\tstart={start}\tstep='
                    f'{1}\tspan={stop-start}\n'
                    )
                for score in scores[:, 1]:
                    out.write(f'{score}\n')

        else:   # unaccepted file type
            raise ValueError(
                'output_file can only have suffixes .wig or .wig.gz.'
                )

    elif (output_file is not None):
        raise TypeError(
            f'output_file is unsupported type "{type(input_file)}". '
            'output_file should be a string specifying the path of the file '
            'to output scores to.'
            )

    if (verbose):
        end_time = time.time()
        print(f'wps took {end_time - start_time} s to complete')

    return scores


def _agg_wps_single_contig(input_file: Union[str, str],
                           contig: str,
                           site_bed: str,
                           window_size: int=120,
                           size_around_sites: int=5000,
                           fraction_low: int=120,
                           fraction_high: int=180,
                           quality_threshold: int=15,
                           verbose: Union[int, bool]=0
                           ):
    """
    Helper function for aggregate_wps. Aggregates wps over sites in one
    contig.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM or SAM file containing paired-end fragment reads or its
        path. `AlignmentFile` must be opened in read mode.
    contig : str
    window_size : int, optional
        Size of window to calculate WPS. Default is k = 120, equivalent
        to L-WPS.
    quality_threshold : int, optional
    workers : int, optional
    verbose : int or bool, optional

    Returns
    -------
    scores : numpy.ndarray
        np array of shape (window_size, 2) where column 1 is the coordinate and
        column 2 is the score.
    """
    if verbose:
        print(f'Aggregating over contig {contig}...')

    # Create tempfile and write contig fragments to
    print(f'Creating frag bed for {contig}')
    _, frag_bed= tf.mkstemp(suffix='.bed.gz', text=True)
    frag_bam_to_bed(input_file,
                    frag_bed,
                    contig=None,
                    quality_threshold=15,
                    verbose=False)

    scores = np.zeros((size_around_sites, 2))

    # Values to add to center of each site to get start and stop of each
    # wps function
    left_of_site = round(-size_around_sites / 2)
    right_of_site = round(size_around_sites / 2)

    assert right_of_site - left_of_site == size_around_sites
    scores[:, 0] = np.arange(left_of_site, right_of_site)

    unaggregated_scores = []

    if (verbose):
        print(f'Opening {input_file} for {contig}...')

    if (verbose >= 2):
        with open(site_bed, 'rt') as sites:
            print('File opened! counting lines for {contig}')
            bed_length = 0
            for line in sites:
                bed_length += 1
    with open(site_bed, 'rt') as sites:
        # verbose stuff
        if (verbose):
            print(f'File opened! Iterating through sites for {contig}...')

        # aggregate wps over sites in bed file
        for line in (
            tqdm(sites, total=bed_length) if verbose>=2 else sites
            ):
            line_items = line.split()
            if ('.' in line_items[5] or contig not in line_items[0]):
                continue
            single_scores = wps(frag_bed,
                                line_items[0],
                                int(line_items[1]) + left_of_site,
                                int(line_items[1]) + right_of_site,
                                output_file=None,
                                window_size=window_size,
                                fraction_low=fraction_low,
                                fraction_high=fraction_high,
                                quality_threshold=quality_threshold,
                                verbose=(verbose-2 if verbose-2>0 else 0)
                                )[:, 1]

            if ('+' in line_items[5]):
                unaggregated_scores.append(single_scores)
            elif ('-' in line_items[5]):
                single_scores = np.flip(single_scores)
                unaggregated_scores.append(single_scores)
            else:   # sites without strand direction are ignored
                pass
        scores[:, 1] = np.sum(unaggregated_scores, axis=0)

        if (verbose):
            print(f'Aggregation complete for {contig}!')

    return scores


def _contig_site_bams(site_bed: str,
                      genome_path: str
                      ) -> dict[str, BinaryIO]:
    with open(genome_path) as genome:
        contigs = [line.split()[0] for line in genome.readlines()]
    tempdir = tf.TemporaryDirectory()
    tempfiles = [f'{tempdir.name}/{contig}.bed' for contig in contigs]
    contig_dict = dict(zip(contigs, tempfiles))
    # TODO: this is known to be broken
    with open(site_bed) as sites:
        for line in sites:
            contig = line.split()[0]
            contig_dict[contig].write(line)
        print(contig_dict.values())
    for file in contig_dict.values():
        file.seek(0)
    return contig_dict


def aggregate_wps(input_file: Union[pysam.AlignmentFile, str],
                  site_bed: str,
                  output_file: str=None,
                  window_size: int=120,
                  size_around_sites: int=5000,
                  fraction_low: int=120,
                  fraction_high: int=180,
                  quality_threshold: int=15,
                  workers: int=1,
                  verbose: Union[bool, int]=0
                  ) -> np.ndarray[int, int]:
    """
    Function that aggregates WPS over sites in BED file
    """
    if (verbose):
        start_time = time.time()
        print(
            f"""
            Calculating aggregate WPS
            input_file: {input_file}
            site_bed: {site_bed}
            output_file: {output_file}
            window_size: {window_size}
            size_around_sites: {size_around_sites}
            quality_threshold: {quality_threshold}
            verbose: {verbose}
            """
            )

    with open(site_bed) as bed_file:
        contigs = []
        for line in bed_file:
            contig = line.split()[0].strip()
            if contig not in contigs:
                contigs.append(contig)

        num_contigs = len(contigs)

    if (verbose):
        print(f'Fragments for {num_contigs} contigs detected.')

    if (verbose >= 2):
        for contig in contigs:
            print(contig)

    input_tuples = zip([input_file] * num_contigs,
                       contigs,
                       [site_bed] * num_contigs,
                       [window_size] * num_contigs,
                       [size_around_sites] * num_contigs,
                       [fraction_low] * num_contigs,
                       [fraction_high] * num_contigs,
                       [quality_threshold] * num_contigs,
                       [verbose - 1 if verbose >= 1 else 0] * num_contigs)

    if (verbose):
        print('Calculating...')

    with Pool(workers) as pool:
        contig_scores = pool.starmap(_agg_wps_single_contig, input_tuples)

    if (verbose):
        print('Compiling scores')

    scores = np.zeros((size_around_sites, 2))
    left_of_site = round(-size_around_sites / 2)
    right_of_site = round(size_around_sites / 2)
    assert right_of_site - left_of_site == size_around_sites
    scores[:, 0] = np.arange(left_of_site, right_of_site)

    for contig_score in contig_scores:
        scores[:, 1] = scores[:, 1] + contig_score[:, 1]

    """
    scores = np.zeros((size_around_sites, 2))

    # Values to add to center of each site to get start and stop of each
    # wps function
    left_of_site = round(-size_around_sites / 2)
    right_of_site = round(size_around_sites / 2)

    assert right_of_site - left_of_site == size_around_sites
    scores[:, 0] = np.arange(left_of_site, right_of_site)

    unaggregated_scores = []

    if (verbose):
        print(f'Opening {input_file}...')

    with pysam.AlignmentFile(input_file) as file:
        if (verbose >= 2):
            with open(site_bed, 'rt') as sites:
                print('File opened! counting lines')
                bed_length = 0
                for line in sites:
                    bed_length += 1
        with open(site_bed, 'rt') as sites:
            # verbose stuff
            if (verbose):
                print(f'File opened! Iterating through sites...')

            # aggregate wps over sites in bed file
            for line in (
                tqdm(sites, total=bed_length) if verbose>=2 else sites
                ):
                line_items = line.split()
                if ('.' in line_items[5]):
                    continue
                single_scores = wps(file,
                                    line_items[0],
                                    int(line_items[1]) + left_of_site,
                                    int(line_items[1]) + right_of_site,
                                    output_file=None,
                                    window_size=window_size,
                                    fraction_low=fraction_low,
                                    fraction_high=fraction_high,
                                    quality_threshold=quality_threshold,
                                    verbose=(verbose-2 if verbose-2>0 else 0)
                                    )[:, 1]

                if ('+' in line_items[5]):
                    unaggregated_scores.append(single_scores)
                elif ('-' in line_items[5]):
                    single_scores = np.flip(single_scores)
                    unaggregated_scores.append(single_scores)
                else:   # sites without strand direction are ignored
                    pass
            scores[:, 1] = np.sum(unaggregated_scores, axis=0)

    if (verbose):
        print(f'Aggregation complete!')
    """


    if (type(output_file) == str):   # check if output specified
        if (verbose):
            print(f'Output file {output_file} specified. Opening...')
        if output_file.endswith(".wig.gz"): # zipped wiggle
            with gzip.open(output_file, 'wt') as out:
                if (verbose):
                    print(f'File opened! Writing...')

                # declaration line
                out.write(
                    f'fixedStep\tchrom=.\tstart={left_of_site}\tstep={1}\tspan'
                    f'={window_size}\n'
                    )
                for score in (tqdm(scores[:, 1])
                              if verbose >= 2
                              else scores[:, 1]):
                    out.write(f'{score}\n')

        elif output_file.endswith(".wig"):  # wiggle
            with open(output_file, 'wt') as out:
                if (verbose):
                    print(f'File opened! Writing...')
                # declaration line
                out.write(
                    f'fixedStep\tchrom=.\tstart={left_of_site}\tstep={1}\tspan'
                    f'={window_size}\n'
                    )
                for score in (tqdm(scores[:, 1])
                              if verbose >= 2
                              else scores[:, 1]):
                    out.write(f'{score}\n')

        else:   # unaccepted file type
            raise ValueError(
                'output_file can only have suffixes .wig or .wig.gz.'
                )

    elif (output_file is not None):
        raise TypeError(
            f'output_file is unsupported type "{type(input_file)}". '
            'output_file should be a string specifying the path of the file '
            'to output scores to.'
            )

    if (verbose):
        end_time = time.time()
        print(f'aggregate_wps took {end_time - start_time} s to complete')

    return scores


# TODO: look through argparse args and fix them all
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculates fragmentation features given a CRAM/BAM/SAM '
        'file',
        epilog='')
    subparsers = parser.add_subparsers(title='subcommands',
                                       dest='subcommand')

    # Common arguments

    # Subcommand 1: frag-coverage
    parser_command1 = subparsers.add_parser('frag-center-coverage',
                                            description=(
                                                'Calculates fragmentation '
                                                'coverage over a region given '
                                                'a CRAM/BAM/SAM file')
                                            )
    # inclusive location of region start in 0-based coordinate system.
    # If not included, will end at the end of the chromosome
    parser_command1.add_argument('--start', type=int)
    # exclusive location of region end in 0-based coordinate system.
    # If not included, will end at the end of the chromosome
    parser_command1.add_argument('--stop', type=int)
    parser_command1.add_argument('--region')   # samtools region string
    parser_command1.add_argument('--method', default="frag-center")
    parser_command1.add_argument('-v', '--verbose', default=False, type=bool)
    parser_command1.set_defaults(func=frag_center_coverage)

    # Subcommand 2: frag-length
    parser_command2 = subparsers.add_parser(
        'frag-length', prog='frag-length',
        description='Calculates fragment lengths given a CRAM/BAM/SAM file'
        )
    parser_command2.add_argument('input_file')
    parser_command2.add_argument('--contig')
    parser_command2.add_argument('--output_file')
    parser_command2.add_argument('--workers', default=1, type=int)
    parser_command2.add_argument('--quality_threshold', default=15, type=int)
    parser_command2.add_argument('-v', '--verbose', action='store_true')
    parser_command2.set_defaults(func=frag_length)

    # Subcommand 3: wps()
    parser_command3 = subparsers.add_parser(
        'wps', prog='wps',
        description='Calculates Windowed Protection Score over a region given '
        'a CRAM/BAM/SAM file'
        )
    parser_command3.add_argument('input_file')
    parser_command3.add_argument('contig')
    # inclusive location of region start in 0-based coordinate system.
    # If not included, will start at the beginning of the chromosome
    parser_command3.add_argument('start', type=int)
    # exclusive location of region end in 0-based coordinate system.
    # If not included, will end at the end of the chromosome
    parser_command3.add_argument('stop', type=int)
    parser_command3.add_argument('-o', '--output_file')
    parser_command3.add_argument('--window_size', default=120, type=int)
    parser_command3.add_argument('-lo', '--fraction_low', default=120,
                                 type=int)
    parser_command3.add_argument('-hi', '--fraction_high', default=180,
                                 type=int)
    parser_command3.add_argument('--quality_threshold', default=15, type=int)
    parser_command3.add_argument('-v', '--verbose', action='count')
    parser_command3.set_defaults(func=wps)

    # Subcommand 4: aggregate-wps
    parser_command4 = subparsers.add_parser(
        'aggregate-wps',
        prog='aggregate-wps',
        description='Calculates Windowed Protection Score over a region '
        'around sites specified in a BED file from alignments in a '
        'CRAM/BAM/SAM file'
        )
    parser_command4.add_argument('input_file')
    parser_command4.add_argument('site_bed')
    parser_command4.add_argument('-o', '--output_file')
    parser_command4.add_argument('--size_around_sites', default=5000, type=int)
    parser_command4.add_argument('--window_size', default=120, type=int)
    parser_command4.add_argument('-lo', '--fraction_low', default=120,
                                 type=int)
    parser_command4.add_argument('-hi', '--fraction_high', default=180,
                                 type=int)
    parser_command4.add_argument('--workers', default=1, type=int)
    parser_command4.add_argument('-v', '--verbose', action='count')
    parser_command4.set_defaults(func=aggregate_wps)


    args = parser.parse_args()
    function = args.func
    funcargs = vars(args)
    funcargs.pop('func')
    funcargs.pop('subcommand')
    # print(funcargs)
    function(**funcargs)