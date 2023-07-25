from __future__ import annotations
import sys
import time
from typing import Union, Tuple

from multiprocessing import Pool

import pysam
import gzip
from numba import jit
from tqdm import tqdm

from finaletools.utils.utils import (
    _not_read1_or_low_quality, _get_contigs, _get_intervals
)


def single_coverage(
        input_file: Union[str, pysam.AlignmentFile],
        contig: str=None,
        start: int=0,
        stop: int=None,
        name: str='.',
        quality_threshold: int=30,
        verbose: Union[bool, int]=False
    ) -> Tuple[str, int, int, str, float]:
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
    name : str
    quality_threshold : int, optional
    verbose : bool, optional

    Returns
    -------
    coverage : int
        Fragment coverage over contig and region.
    """
    # TODO: determine if reference (like as found in pysam) is necessary
    # TODO: consider including region string (like in pysam)
    if verbose:
        start_time = time.time()

    # initializing variable for coverage tuple outside of with statement
    coverage = 0

    input_type = type(input_file)

    try:
        # Handling different input types
        if input_type == pysam.AlignmentFile:
            sam_file = input_file
        elif input_type == str:
            sam_file = pysam.AlignmentFile(input_file, 'r')
        else:
            raise TypeError('input_file should be str or pysam.AlignmentFile')

        # Iterating on each read in file in
        # specified contig/chromosome
        for read1 in sam_file.fetch(
            contig=contig,
            start=0 if (tempstart:=start-500) < 0 else tempstart,
            stop=stop+500 if stop is not None else None
        ):
            # Only select forward strand and filter out
            # non-paired-end reads and low-quality reads
            if _not_read1_or_low_quality(read1, quality_threshold):
                pass
            else:
                # calculate mid-point of fragment
                center = read1.reference_start + read1.template_length // 2
                if start is None and stop is None:
                    coverage += 1
                elif start is None:
                    if center < stop:
                        coverage += 1
                elif stop is None:
                    if center >= start:
                        coverage += 1
                elif (center >= start) and (center < stop):
                    coverage += 1
                else:
                    pass
    finally:
        if input_type == str:
            sam_file.close()


    if verbose:
        end_time = time.time()
        tqdm.write(
            f'frag_coverage took {end_time - start_time} s to complete\n'
        )

    return contig, start, stop, name, coverage


def _single_coverage_star(args):
    """
    Helper function that takes a tuple of args and applies to
    single_coverage. To be used in imap.
    """
    return single_coverage(*args)


def coverage(
        input_file: Union[str, pysam.AlignmentFile],
        interval_file: str,
        output_file: str,
        scale_factor: float=1e6,
        quality_threshold: int=30,
        workers: int=1,
        verbose: Union[bool, int]=False):
    """
    Return estimated fragment coverage over intervals specified in
    `intervals`. Fragments are read from `input_file` which may be
    either a BAM or SAM file. Uses an algorithm where the midpoints of
    fragments are calculated and coverage is tabulated from the
    midpoints that fall into the specified region. Not suitable for
    fragments of size approaching interval size.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, SAM, or CRAM file containing paired-end fragment reads or
        its path. `AlignmentFile` must be opened in read mode.
    intervals : str
        Path for BAM file containing intervals
    output_file : string, optional
        Path for bed file to print coverages to. If output_file = `_`,
        results will be printed to stdout.
    scale_factor : int, optional
        Amount to multiply coverages by. Default is 10^6.
    quality_threshold : int, optional
    verbose : int or bool, optional

    Returns
    -------
    coverage : int
        Fragment coverage over contig and region.
    """
    if (verbose):
        start_time = time.time()
        sys.stderr.write(
            f"""
            input_file: {input_file}
            interval file: {interval_file}
            output_file: {output_file}
            quality_threshold: {quality_threshold}
            workers: {workers}
            verbose: {verbose}\n
            """
        )

    contigs = _get_contigs(
        input_file,
        verbose=verbose-1 if verbose>1 else 0
    )

    # calculating coverage for each contig in the input_file
    contig_intervals = []
    for contig, length in contigs:
        contig_intervals.append((
            input_file,
            contig,
            0,
            length,
            ".",
            quality_threshold,
            verbose - 1 if verbose > 1 else 0
        ))
    if verbose:
        sys.stderr.write('Creating process pool\n')
    try:
        pool = Pool(processes=workers, )
        if verbose:
            sys.stderr.write('Calculating total coverage for file,\n')

        total_coverage_results = pool.imap(
            _single_coverage_star,
            tqdm(
                contig_intervals,
                desc='Genome contigs',
                position=0
            ) if verbose else contig_intervals
        )

        if verbose:
            tqdm.write('reading intervals\n')

        intervals = pool.apply(
            _get_intervals,
            (input_file, interval_file, quality_threshold, verbose)
        )

        if verbose:
            tqdm.write('calculating coverage\n')

        coverages = pool.imap(
            _single_coverage_star,
            tqdm(
                intervals,
                desc='Intervals',
                position=2
            ) if verbose else intervals,
            len(intervals) // 2 // workers + 1
        )

        if verbose:
            tqdm.write('Retrieving total coverage for file\n')
        total_coverage = sum(coverage[4] for coverage in total_coverage_results)
        if verbose:
                tqdm.write(f'Total coverage is {total_coverage}\n')

        # Output
        output_is_file = False

        if output_file != None:
            if verbose:
                tqdm.write('Writing results to output\n')
            try:
                # handle output types
                if output_file.endswith('.bed'):
                    output_is_file = True
                    output = open(output_file, 'w')
                elif output_file.endswith('.bed.gz'):
                    output = gzip.open(output_file, 'w')
                    output_is_file = True
                elif output_file == '-':
                    output = sys.stdout
                else:
                    raise ValueError(
                        'output_file should have .bed or .bed.gz as as suffix'
                    )

                # print to files
                for coverage in coverages:
                    output.write(
                        f'{coverage[0]}\t{coverage[1]}\t{coverage[2]}\t'
                        f'{coverage[3]}\t'
                        f'{coverage[4]/total_coverage*scale_factor}\n'
                    )

            finally:
                if output_is_file:
                    output.close()
    finally:
        pool.close()

    if verbose:
        end_time = time.time()
        sys.stderr.write(
            f'coverage took {end_time - start_time} s to complete\n'
        )

    return coverages

