from __future__ import annotations
import sys
import time
from typing import Union
from pathlib import Path
from multiprocessing import Pool
from functools import partial

import pysam
import gzip
from tqdm import tqdm

from finaletoolkit.utils.utils import (
    _get_intervals, frag_generator
)


def single_coverage(
    input_file: Union[str, pysam.TabixFile, pysam.AlignmentFile, Path],
    contig: str | None = None,
    start: int | None = 0,
    stop: int | None = None,
    name: str | None = '.',
    min_length: int | None = None,
    max_length: int | None = None,
    intersect_policy: str = "midpoint",
    quality_threshold: int = 30,
    verbose: bool | int = False
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
        Name of interval. Default is '.'. '.' is used if `None` is
        supplied.
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
        min_length=min_length,
        max_length=max_length,
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
    
    adjusted_name = '.' if name is None else name

    return contig, start, stop, adjusted_name, coverage


def _single_coverage_star(partial_coverage, interval):
    contig, start, stop, name = interval
    return partial_coverage(contig=contig, start=start, stop=stop, name=name)


def coverage(
        input_file: Union[str, pysam.TabixFile, pysam.AlignmentFile, Path],
        interval_file: str,
        output_file: str,
        scale_factor: float=1.,
        min_length: int | None=None,
        max_length: int | None=None,
        normalize: bool=False,
        intersect_policy: str="midpoint",
        quality_threshold: int=30,
        workers: int=1,
        verbose: Union[bool, int]=False
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
    coverages : Iterable[tuple[str, int, int, str, float]]
        Fragment coverages over intervals.
    """
    returnVal = []
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
        pool = Pool(processes=workers, )
        if verbose:
            tqdm.write('Calculating total coverage for file,\n')

        if verbose:
            tqdm.write('reading intervals\n')

        if normalize:
            # queue in pool
            total_coverage_results = pool.apply_async(
                single_coverage,
                (input_file, None, 0, None, '.'),
                {"min_length": min_length,
                 "max_length": max_length,
                 "intersect_policy": intersect_policy,
                 "quality_threshold": quality_threshold,
                 "verbose": verbose}
            )

        intervals = _get_intervals(
            interval_file)

        if verbose:
            tqdm.write('calculating coverage\n')

        partial_single_coverage = partial(
            single_coverage, input_file=input_file, min_length=min_length,
            max_length=max_length, intersect_policy=intersect_policy,
            quality_threshold=quality_threshold, verbose=max(0, verbose-1))
        coverages = pool.imap(
            partial(_single_coverage_star, partial_single_coverage), intervals,
            chunksize=max(len(intervals)//workers, 1))

        if verbose:
            tqdm.write('Retrieving total coverage for file\n')

        # normalize
        if normalize:
            total_coverage = total_coverage_results.get()
            if verbose:
                    tqdm.write(f'Total coverage is {total_coverage}\n')
            scale_factor /= total_coverage[4]

        # Output
        output_is_file = False

        if output_file != None:
            if verbose:
                tqdm.write('Writing results to output\n')
            try:
                # handle output types
                if (output_file.endswith('.bed')
                    or output_file.endswith('.bedgraph')):
                    output_is_file = True
                    output = open(output_file, 'w')
                elif output_file.endswith('.bed.gz'):
                    output = gzip.open(output_file, 'wt')  # text writing
                    output_is_file = True
                elif output_file == '-':
                    output = sys.stdout
                else:
                    raise ValueError(
                        'output_file should have .bed or .bed.gz as suffix')
            
                # print to files
                if output_file.endswith(".bedgraph"):
                    for contig, start, stop, name, coverage in coverages:
                        output.write(
                            f'{contig}\t{start}\t{stop}\t'
                            f'{coverage * scale_factor}\n'
                        )
                        returnVal.append(
                            (contig, start, stop, name,
                             coverage * scale_factor))
                else:
                    for contig, start, stop, name, coverage in coverages:
                        output.write(
                            f'{contig}\t{start}\t{stop}\t'
                            f'{name}\t'
                            f'{coverage * scale_factor}\n'
                        )
                        returnVal.append(
                            (contig, start, stop, name,
                             coverage * scale_factor))
            finally:
                if output_is_file:
                    output.close()
        else:
            returnVal = [(contig, start, stop, name, coverage * scale_factor)
                         for (contig, start, stop, name, coverage)
                         in coverages]
    finally:
        pool.close()

    if verbose:
        end_time = time.time()
        tqdm.write(
            f'coverage took {end_time - start_time} s to complete\n'
        )
    return returnVal

