from __future__ import annotations
import sys
import gzip
import time
import pysam
from tqdm import tqdm
from typing import Union, Tuple
from multiprocessing import Pool
from functools import partial
from finaletools.utils.utils import _get_intervals, frag_generator

def single_coverage(
        input_file: Union[str, pysam.AlignmentFile],
        contig: str,
        start: int,
        stop: int,
        name: str,
        min_length: int,
        max_length: int,
        intersect_policy: str,
        quality_threshold: int,
        verbose: bool
    ) -> Tuple[str, int, int, str, float]:

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

    frags = frag_generator(
        input_file=input_file,
        contig=contig,
        quality_threshold=quality_threshold,
        start=start,
        stop=stop,
        min=min_length,
        max=max_length,
        intersect_policy=intersect_policy,
        verbose=verbose)

    coverage = sum(1 for _ in frags)
        
    if verbose:
        end_time = time.time()
        tqdm.write(
            f'Calculating coverage took {end_time - start_time} s to complete.\n'
        )

    return contig, start, stop, name, coverage

def _single_coverage_star(partial_coverage, interval):
    contig, start, stop, name = interval
    return partial_coverage(contig=contig, start=start, stop=stop, name=name)

def coverage_intervals(
        input_file: Union[str, pysam.AlignmentFile],
        interval_file: str,
        output_file: str,
        scale_factor: float,
        min_length: int,
        max_length: int,
        intersect_policy: str,
        quality_threshold: int,
        workers: int,
        normalize: bool,
        verbose: bool
    ):

    if (verbose):
        start_time = time.time()
        sys.stderr.write(
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
        sys.stderr.write('Creating process pool.\n')
    try:
        pool = Pool(processes=workers)
        if verbose:
            sys.stderr.write('Calculating total coverage for file.\n')

        if normalize:
            total_coverage_results = pool.apply_async(single_coverage, (input_file, None, 0, None, '.', min_length, max_length, intersect_policy, quality_threshold, False))

        if verbose:
            tqdm.write('Calculating coverage...\n')

        intervals = _get_intervals(interval_file)

        partial_single_coverage = partial(single_coverage, input_file=input_file, min_length=min_length, max_length=max_length, intersect_policy=intersect_policy, quality_threshold=quality_threshold, verbose=verbose)
        coverages = pool.imap(partial(_single_coverage_star, partial_single_coverage), intervals, chunksize=max(len(intervals)//workers, 1))

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

        if output_file != None:
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
                    for contig, start, stop, name, coverage in coverages:
                        output.write(
                            f'{contig}\t{start}\t{stop}\t'
                            f'{coverage/total_coverage*scale_factor}\n'
                        )
                else:
                    for contig, start, stop, name, coverage in coverages:
                        output.write(
                            f'{contig}\t{start}\t{stop}\t'
                            f'{name}\t'
                            f'{coverage/total_coverage*scale_factor}\n'
                        )

            finally:
                if output_is_file:
                    output.close()
    finally:
        pool.close()

    if verbose:
        end_time = time.time()
        sys.stderr.write(
            f'Coverage over intervals took {end_time - start_time} s to complete\n'
        )

    return coverages

