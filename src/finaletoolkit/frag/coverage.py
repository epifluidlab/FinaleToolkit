from __future__ import annotations
import sys
import time
from typing import Union, Tuple

from multiprocessing import Pool

import pysam
import gzip
from tqdm import tqdm

from finaletoolkit.utils.utils import (
    _not_read1_or_low_quality, _get_contigs, _get_intervals, frag_generator
)


def single_coverage(
        input_file: Union[str, pysam.AlignmentFile],
        contig: str=None,
        start: int=0,
        stop: int=None,
        name: str='.',
        intersect_policy: str="midpoint",
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
        BAM, SAM, CRAM, or Frag.gz file containing paired-end fragment reads or
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
        tqdm.write(
            f"""
input_file: 
contig: {contig}
start: {start}
stop: {stop}
name: {name}
intersect_policy: {intersect_policy}
quality_threshold: {quality_threshold}
verbose: {verbose}
"""
        )

    # initializing variable for coverage tuple outside of with statement
    coverage = 0

    input_type = type(input_file)

    frags = frag_generator(
        input_file=input_file,
        contig=contig,
        quality_threshold=quality_threshold,
        start=start,
        stop=stop,
        fraction_low=10,
        fraction_high=10000000,
        intersect_policy=intersect_policy)

    # Iterating on each frag in file in
    # specified contig/chromosome
    for frag in frags:
        coverage += 1
        
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
        verbose: Union[bool, int]=False
    ):
    """
    Return estimated fragment coverage over intervals specified in
    `intervals`. Fragments are read from `input_file` which may be
    a SAM, BAM, CRAM, or Frag.gz file. Uses an algorithm where the midpoints of
    fragments are calculated and coverage is tabulated from the
    midpoints that fall into the specified region. Not suitable for
    fragments of size approaching interval size.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        SAM, BAM, CRAM, or Frag.gz file containing paired-end fragment 
        reads or its path. `AlignmentFile` must be opened in read mode.
    interval_file : str
        BED4 file containing intervals over which to generate coverage
        statistics.
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
    #FIXME update docstring
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

    if verbose:
        sys.stderr.write('Creating process pool\n')
    try:
        pool = Pool(processes=workers, )
        if verbose:
            sys.stderr.write('Calculating total coverage for file,\n')

        total_coverage_results = pool.apply_async(
            single_coverage,
            (input_file, None, 0, None, '.', "midpoint", quality_threshold, False)
        )

        if verbose:
            tqdm.write('reading intervals\n')

        intervals = _get_intervals(
            input_file, interval_file,
            intersect_policy="midpoint",
            quality_threshold=quality_threshold,
            verbose=verbose)

        if verbose:
            tqdm.write('calculating coverage\n')

        coverages = pool.imap(
            _single_coverage_star,
            tqdm(
                intervals,
                desc='Intervals',
                position=2
            ) if verbose else intervals,
            min(len(intervals) // 2 // workers + 1, 40)
        )

        if verbose:
            tqdm.write('Retrieving total coverage for file\n')
        total_coverage = total_coverage_results.get()
        if verbose:
                tqdm.write(f'Total coverage is {total_coverage}\n')

        # Output
        output_is_file = False

        if output_file != None:
            if verbose:
                tqdm.write('Writing results to output\n')
            try:
                # handle output types
                if (output_file.endswith('.bed')
                    or output_file.endswith('.bedgraph')
                ):
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
                if output_file.endswith(".bedgraph"):
                    for contig, start, stop, name, coverage in coverages:
                        output.write(
                            f'{contig}\t{start}\t{stop}\t'
                            f'{coverage/total_coverage[4]*scale_factor}\n'
                        )
                else:
                    for contig, start, stop, name, coverage in coverages:
                        output.write(
                            f'{contig}\t{start}\t{stop}\t'
                            f'{name}\t'
                            f'{coverage/total_coverage[4]*scale_factor}\n'
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

