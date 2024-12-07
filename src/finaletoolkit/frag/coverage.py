from __future__ import annotations
import sys
import time
from typing import Union, Iterable
from pathlib import Path
from multiprocessing import Pool
from functools import partial

import pysam
import gzip
from tqdm import tqdm

from finaletoolkit.utils.utils import (
    _get_intervals, frag_generator
)
from finaletoolkit.utils._deprecation import deprecated


def _single_coverage(
        input_file: Union[str, pysam.TabixFile, pysam.AlignmentFile, Path],
        contig: str | None,
        start: int | None,
        stop: int | None,
        name: str,
        min_length: int | None,
        max_length: int | None,
        scale_factor: float,
        intersect_policy: str,
        quality_threshold: int,
        verbose: Union[bool, int]
    ) -> tuple[str | None, int | None, int | None, str, float]:
    """
    Return estimated fragment coverage over specified `contig` and
    region of`input_file`. Uses an algorithm where the midpoints of
    fragments are calculated and coverage is tabulated from the
    midpoints that fall into the specified region. Not suitable for
    fragments of size approaching region size.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, SAM, CRAM, or Frag.gz file containing paired-end fragment
        reads or its path. `AlignmentFile` must be opened in read mode.
    contig : string, optional
        Default is None.
    start : int, optional
        0-based start coordinate to get coverage from. Default is 0.
    stop : int, optional
        1-based stop coordinate to get coverage from. Default is None
        and goes to end of chromosome/contig.
    name : str, optional
        Name of interval. Default is '.'.
    min_length: int or None, optional
        Minimum length of fragments to be included.
    max_length: int or None, optional
        Maximum length of fragments to be included.
    intersect_policy: str, optional
        Specifies how to determine whether fragments are in interval.
        'midpoint' (default) calculates the central coordinate of each
        fragment and only selects the fragment if the midpoint is in the
        interval. 'any' includes fragments with any overlap with the
        interval.
    quality_threshold : int, optional
        Minimum MAPQ to filter for. Default is 30.
    verbose : bool, optional
        Prints messages to stderr. Default is false.
    Returns
    -------
    (contig, start, stop, name, coverage) : tuple[str, int, int, str, float]
        Fragment coverage over contig and region.
    """
    if verbose:
        start_time = time.time()
        tqdm.write(
            f"""
            input_file: {input_file}
            contig: {contig}
            start: {start}
            stop: {stop}
            name: {name}
            min_length: {min_length}
            max_length: {max_length}
            intersect_policy: {intersect_policy}
            quality_threshold: {quality_threshold}
            verbose: {verbose}
            """
        )

    # initializing variable for coverage tuple outside of with statement
    coverage = 0

    frags = frag_generator(
        input_file=input_file,
        contig=contig,
        quality_threshold=quality_threshold,
        start=start,
        stop=stop,
        fraction_low=min_length,
        fraction_high=max_length,
        intersect_policy=intersect_policy)

    # Iterating on each frag in file in
    # specified contig/chromosome
    for frag in frags:
        coverage += 1
        
    if verbose:
        end_time = time.time()
        tqdm.write(
            f'single_coverage took {end_time - start_time} s to complete\n'
        )

    return contig, start, stop, name, coverage*scale_factor


def _single_coverage_star(partial_coverage, interval):
    contig, start, stop, name = interval
    return partial_coverage(contig=contig, start=start, stop=stop, name=name)


def coverage(
        input_file: Union[str, pysam.TabixFile, pysam.AlignmentFile, Path],
        interval_file: str,
        output_file: str=None,
        scale_factor: float=1.,
        min_length: int | None=None,
        max_length: int | None=None,
        normalize: bool=False,
        intersect_policy: str="midpoint",
        quality_threshold: int=30,
        workers: int=1,
        verbose: Union[bool, int]=False,
    ) -> list[tuple[str, int, int, str, float]]:
    """
    Return estimated fragment coverage over intervals specified in
    `intervals`. Fragments are read from `input_file` which may be
    a SAM, BAM, CRAM, or Frag.gz file. Uses an algorithm where the
    midpoints of fragments are calculated and coverage is tabulated from
    the midpoints that fall into the specified region. Not suitable for
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
        Path for bed file to print coverages to. If output_file = `-`,
        results will be printed to stdout.
    scale_factor : int, optional
        Amount to multiply coverages by. Default is 1.
    min_length: int or None, optional
        Minimum length of fragments to be included.
    max_length: int or None, optional
        Maximum length of fragments to be included.
    normalize : bool
        When set to True, divide by total coverage
    intersect_policy: str, optional
        Specifies how to determine whether fragments are in interval.
        'midpoint' (default) calculates the central coordinate of each
        fragment and only selects the fragment if the midpoint is in the
        interval. 'any' includes fragments with any overlap with the
        interval.
    quality_threshold : int, optional
        Minimum MAPQ. Default is 30.
    workers : int, optional
        Number of subprocesses to spawn. Increases speed at the expense
        of memory.
    verbose : int or bool, optional

    Returns
    -------
    coverages : Iterable or list of (str, int, int, str, float)
        Fragment coverages over intervals.
    """
    if (verbose):
        start_time = time.time()
        tqdm.write(
            f"""
            input_file: {input_file}
            interval file: {interval_file}
            output_file: {output_file}
            scale_factor: {scale_factor}
            min_length: {min_length}
            max_length: {max_length}
            intersect_policy: {intersect_policy}
            quality_threshold: {quality_threshold}
            workers: {workers}
            normalize: {normalize}
            verbose: {verbose} \n
            """
        )

    if verbose:
        tqdm.write('Creating process pool\n')
    try:
        pool = Pool(processes=workers)
        if verbose:
            tqdm.write('Calculating total coverage for file,\n')
        # XXX: implementation in develop branch
        """
        if normalize:
            # queue in pool
            total_coverage_results = pool.apply_async(
                single_coverage,
                (input_file, None, 0, None, '.'),
                {"intersect_policy": "midpoint",
                 "quality_threshold": quality_threshold,
                 "verbose": False}
            )

        intervals = _get_intervals(interval_file)
        """

        if normalize:
            total_coverage_results = pool.apply_async(
                _single_coverage,
                (input_file, None, 0, None, '.', min_length, max_length,
                 1., intersect_policy, quality_threshold, False))

        if verbose:
            tqdm.write('Calculating coverage...\n')

        intervals = _get_intervals(interval_file)

        partial_single_coverage = partial(
            _single_coverage, input_file=input_file, min_length=min_length,
            scale_factor=scale_factor, max_length=max_length,
            intersect_policy=intersect_policy,
            quality_threshold=quality_threshold, verbose=verbose)
        coverages = pool.imap(
            partial(_single_coverage_star, partial_single_coverage),
            intervals, chunksize=max(len(intervals)//workers, 1))
        
        if verbose:
            tqdm.write('Retrieving total coverage for file\n')
        
        if normalize:
            total_coverage = total_coverage_results.get()
            total_coverage = total_coverage[4]
        else:
            total_coverage = 1

        if verbose and normalize:
                tqdm.write(f'Total coverage is {total_coverage}\n')
            
        output_is_file = False
        return_val = []
        if output_file is not None:
            coverages_list = [i for i in coverages]
            if verbose:
                tqdm.write('Writing results to output. \n')
            try:
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
                        'The output file should have .bed or .bed.gz as as suffix.'
                    )

                if output_file.endswith(".bedgraph"):
                    for contig, start, stop, name, coverage in coverages_list:
                        output.write(
                            f'{contig}\t{start}\t{stop}\t'
                            f'{coverage * scale_factor}\n'
                        )
                        return_val.append(
                            (contig, start, stop, name,
                             coverage * scale_factor))
                else:
                    for contig, start, stop, name, coverage in coverages_list:
                        output.write(
                            f'{contig}\t{start}\t{stop}\t'
                            f'{name}\t'
                            f'{coverage * scale_factor}\n'
                        )
                        return_val.append(
                            (contig, start, stop, name,
                             coverage * scale_factor))
            finally:
                if output_is_file:
                    output.close()
        else:
            return_val = [(contig, start, stop, name, coverage * scale_factor)
                         for (contig, start, stop, name, coverage)
                         in coverages]
    finally:
        pool.close()

    if verbose:
        end_time = time.time()
        sys.stderr.write(
            f'Coverage over intervals took {end_time - start_time} s to complete\n'
        )
    return return_val

# deprecated functions
@deprecated
def single_coverage(
    input_file: Union[str, pysam.TabixFile, pysam.AlignmentFile, Path],
    contig: str | None=None,
    start: int | None=None,
    stop: int | None=None,
    name: str='.',
    min_length: int | None=None,
    max_length: int | None=None,
    scale_factor=1.,
    intersect_policy: str="midpoint",
    quality_threshold: int=30,
    verbose: Union[bool, int]=False
):
    return _single_coverage(
        input_file=input_file,
        contig=contig,
        start=start,
        stop=stop,
        name=name,
        min_length=min_length,
        max_length=max_length,
        intersect_policy=intersect_policy,
        scale_factor=scale_factor,
        quality_threshold=quality_threshold,
        verbose=verbose,
    )