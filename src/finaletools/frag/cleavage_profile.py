from __future__ import annotations
from typing import Union
from sys import stderr, stdout, stdin
from multiprocessing import Pool
import gzip
import time
import numpy as np
import pyBigWig as pbw
import pysam

from finaletools.utils.utils import frag_array, _parse_chrom_sizes, _reduce_overlaps_in_file, _convert_to_list, _merge_all_intervals

def cleavage_profile(
    input_file: str,
    contig: str,
    start: int,
    stop: int,
    min_length: int,
    max_length: int,
    quality_threshold: int,
    five: bool,
    verbose: bool
) -> np.ndarray:
    if (verbose):
        start_time = time.time()
        stderr.write(
            f"""
            Calculating cleavage profile
            input_file: {input_file}
            contig: {contig}
            start: {start}
            stop: {stop}
            min_length: {min_length}
            max_length: {max_length}
            quality_threshold: {quality_threshold}
            five: {five}
            verbose: {verbose}
            """
        )

    frags = frag_array(
        input_file=input_file,
        contig=contig,
        quality_threshold=quality_threshold,
        start=start,
        stop=stop,
        min=min_length,
        max=max_length,
        intersect_policy="any",
        verbose=verbose
    )

    positions = np.arange(start, stop)
    if (verbose):
        stderr.write('Beginning cleavage profile calculation...\n')
    fragwise_overlaps = np.logical_and(
        np.greater_equal(positions[np.newaxis], frags['start'][:,np.newaxis]),
        np.less(positions[np.newaxis], frags['stop'][:,np.newaxis])
    )
    depth = np.sum(fragwise_overlaps, axis=0)
    if np.sum(depth) == 0:
        return None
    if five:
        forward_ends = np.logical_and(
            np.equal(
                positions[np.newaxis], frags['start'][:, np.newaxis]
            ), frags['strand'][:, np.newaxis]
        )
        reverse_ends = np.logical_and(
            np.equal(
                positions[np.newaxis], frags['stop'][:, np.newaxis]
            ), np.logical_not(frags['strand'][:, np.newaxis])
        )
    elif not five:
        forward_ends = np.logical_and(
            np.equal(
                positions[np.newaxis], frags['stop'][:, np.newaxis]
            ), frags['strand'][:, np.newaxis]
        )
        reverse_ends = np.logical_and(
            np.equal(
                positions[np.newaxis], frags['start'][:, np.newaxis]
            ), np.logical_not(frags['strand'][:, np.newaxis])
        )
    ends = np.sum(np.logical_or(forward_ends, reverse_ends), axis=0)
    proportions = np.zeros_like(depth, dtype=np.float64)
    non_zero_mask = depth != 0
    proportions[non_zero_mask] = ends[non_zero_mask] / depth[non_zero_mask] * 100

    results = np.zeros_like(proportions, dtype=[
        ('contig', 'U16'),
        ('pos', 'i8'),
        ('proportion', 'f8'),
    ])
    results['contig'] = contig
    results['pos'] = np.arange(start, stop)
    results['proportion'] = proportions

    return results


def _cleavage_profile_star(args):
    return cleavage_profile(*args)


def cleavage_profile_intervals(
    input_file: str,
    interval_file: str,
    min_length: int,
    max_length: int,
    quality_threshold: int,
    output_file: str,
    chrom_sizes: str,
    five: bool,
    workers: int,
    verbose: bool
   
):
    if (verbose):
        start_time = time.time()
        stderr.write(
            f"""
            Calculating cleavage profile:
            input_file: {input_file}
            interval_file: {interval_file}
            min_length: {min_length}
            max_length: {max_length}
            quality_threshold: {quality_threshold}
            output_file: {output_file}
            chrom_sizes: {chrom_sizes}
            five: {five}
            workers: {workers}
            verbose: {verbose}
            """
        )
    
    if (input_file == '-' and interval_file == '-'):
        raise ValueError('The input file and interval file cannot both read from stdin')
    
    if (input_file.endswith('.sam')
        or input_file.endswith('.bam')
        or input_file.endswith('.cram')):
        with pysam.AlignmentFile(input_file, 'r') as bam:
            references = bam.references
            lengths = bam.lengths
            header = list(zip(references, lengths))
    elif (input_file.endswith('.bed')
          or input_file.endswith('.bed.gz')
          or input_file.endswith('.frag')
          or input_file.endswith('.frag.gz')
    ):
        if chrom_sizes is None:
            raise ValueError(
                '--chrom_sizes must be specified for BED/Fragment files'
            )
        header = _parse_chrom_sizes(chrom_sizes)
    else:
        raise ValueError("Not a supported file type.")

    if (verbose > 1):
        stderr.write(f'Header is {header}\n')
    
    reduced_intervals = _reduce_overlaps_in_file(interval_file)
    converted_intervals = _convert_to_list(reduced_intervals)
    all_intervals = _merge_all_intervals(converted_intervals)
    contigs, starts, stops = zip(*all_intervals)
    count = len(contigs)


    if (verbose):
        stderr.write('Zipping inputs. \n')

    interval_list = zip(
        count*[input_file],
        contigs,
        starts,
        stops,
        count*[min_length],
        count*[max_length],
        count*[quality_threshold],
        count*[five],
        count*[verbose]
    )

    if (verbose):
        stderr.write('Calculating cleavage profile...\n')

    try:
        pool = Pool(processes=workers)
        interval_scores = pool.imap(_cleavage_profile_star, interval_list, chunksize=max(len(contigs)//workers, 1))

        # output
        if (type(output_file) == str):   # check if output specified
            if (verbose):
                stderr.write(
                    f'Output file {output_file} specified. Opening...\n'
                )

            if output_file.endswith(".bw"):
                with pbw.open(output_file, 'w') as bigwig:
                    bigwig.addHeader(header)
                    for interval_score in interval_scores:
                        if interval_score is not None:
                            contigs = interval_score['contig']
                            starts = interval_score['pos']
                            scores = interval_score['proportion']
                            stops = starts + 1

                            if contigs.shape == (0,):
                                continue

                            try:
                                bigwig.addEntries(
                                    chroms=contigs,
                                    starts=starts,
                                    ends=stops,
                                    values=scores.astype(np.float64),
                                )
                            except RuntimeError as e:
                                stderr.write(
                                    f'{contigs[0]}:{starts[0]}-{stops[-1]}\n'
                                )
                                stderr.write(
                                    '/n Invalid or out of order interval '
                                    'encountered. Skipping to next.\n'
                                )
                                continue
            else:
                raise ValueError(
                    'Output file can only have suffix .bw'
                    )

        elif (output_file is not None):
            raise TypeError(
                f'The output file is unsupported type "{type(output_file)}". '
                'Output file should be a string specifying the path of the .bw'
                'file to output scores to.'
                )
    finally:
        pool.close()

    if (verbose):
        end_time = time.time()
        stderr.write(
            f'Calculating cleavage profile took {end_time - start_time} s to complete\n'
        )
    
    return None
