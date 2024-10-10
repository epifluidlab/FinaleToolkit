"""Cleavage Profiler

This module is used to find the cleavage profile as described in Zhou
et al 2022 (https://doi.org/10.1073/pnas.2209852119). Cleavage profile
describes the proportion of fragment ends at a site over the depth at
the site (cleavage proportion) calculated over a 5+/- nt window around a
CpG site.
"""

from __future__ import annotations
from typing import Union
from pathlib import Path
from sys import stderr
from multiprocessing import Pool
import gzip
import time

import numpy as np
import pyBigWig as pbw
import pysam

from finaletoolkit.utils.utils import (
    frag_array, chrom_sizes_to_list, _reduce_overlaps_in_file,
    _convert_to_list, _merge_all_intervals, chrom_sizes_to_dict
    )


def cleavage_profile(
    input_file: str,
    chrom_size: int,
    contig: str,
    start: int,
    stop: int,
    left: int=0,
    right: int=0,
    fraction_low: int=1,
    fraction_high: int=10000000,
    quality_threshold: int=30,
    verbose: Union[bool, int]=0
) -> np.ndarray:
    """
    Cleavage profile calculated over a single interval.

    Parameters
    ---------
    input_file: str
        SAM, BAM, CRAM, or FRAG file with fragment information.
    chrom_size: int
        length of contig.
    contig: str
        Chromosome or contig
    start: int
        0-based start coordinate
    stop: int
        1-based end coordinate
    left: int
        Amount to subtract from start coordinate. Useful if only given
        coordinates of CpG.
    right: int
        Amount to add to stop coordinate.
    fraction_low: int
        Minimum fragment size to include
    fraction_high: int
        Maximum fragment size to include
    quality_threshold: int
        Minimum MAPQ
    verbose: bool or in

    Return
    ------
    cleavage_proportions: NDArray
        Array of cleavage proportions over given interval.
    """
    if (verbose):
        start_time = time.time()
        stderr.write(
            f"""
            Calculating cleavage profile
            input_file: {input_file}
            contig: {contig}
            start: {start}
            stop: {stop}
            fraction_low: {fraction_low}
            fraction_high: {fraction_high}
            quality_threshold: {quality_threshold}
            verbose: {verbose}
            """
        )
    adj_start = max(start-left, 0)
    adj_stop = min(stop+right, chrom_size)

    frags = frag_array(
        input_file=input_file,
        contig=contig,
        quality_threshold=quality_threshold,
        start=adj_start,
        stop=adj_stop,
        fraction_low=fraction_low,
        fraction_high=fraction_high,
        intersect_policy="any"
    )

    positions = np.arange(adj_start, adj_stop)

    # finding depth at sites
    fragwise_overlaps = np.logical_and(
        np.greater_equal(positions[np.newaxis], frags['start'][:,np.newaxis]),
        np.less(positions[np.newaxis], frags['stop'][:,np.newaxis])
    )
    depth = np.sum(fragwise_overlaps, axis=0)

    # finding ends
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
    results['pos'] = positions
    results['proportion'] = proportions


    return results


def _cleavage_profile_star(args):
    return cleavage_profile(*args)


def _cli_cleavage_profile(
    input_file: Union[str, Path],
    interval_file: Union[str, Path],
    chrom_sizes: Union[str, Path],
    left: int=0,
    right: int=0,
    fraction_low: int=1,
    fraction_high: int=10000000,
    quality_threshold: int=30,
    output_file: str='-',
    workers: int=1,
    verbose: Union[bool, int]=0
):
    """
    Function called when running cleavage profile subcommand in cli.
    Multithreaded implementation over intervals in bed.
    """

    if (verbose):
        start_time = time.time()
        stderr.write(
            f"""
            Calculating cleavage profile
            input_file: {input_file}
            interval_file: {interval_file}
            fraction_low: {fraction_low}
            fraction_high: {fraction_high}
            quality_threshold: {quality_threshold}
            output_file: {output_file}
            workers: {workers}
            verbose: {verbose}
            """
        )
    
    if (input_file == '-' and interval_file == '-'):
        raise ValueError('input_file and site_bed cannot both read from stdin')
    
    if chrom_sizes is None:
        raise ValueError(
            '--chrom_sizes must be specified'
        )
    
      # get header from input_file
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
        header = chrom_sizes_to_list(chrom_sizes)
    else:
        raise ValueError("Not a supported file type.")

    if (verbose > 1):
        stderr.write(f'header is {header}\n')
    
    # reading intervals from bed and removing overlaps
    # NOTE: assumes that bed file is sorted.
    reduced_intervals = _reduce_overlaps_in_file(interval_file)
    converted_intervals = _convert_to_list(reduced_intervals)
    all_intervals = _merge_all_intervals(converted_intervals)
    contigs, starts, stops = zip(*all_intervals)

    # reading chrom.sizes file

    size_dict = chrom_sizes_to_dict(chrom_sizes)

    sizes = [size_dict[contig] for contig in contigs]
    
    count = len(contigs)

    if (verbose):
        stderr.write('Zipping inputs\n')

    interval_list = zip(
        count*[input_file],
        sizes,
        contigs,
        starts,
        stops,
        count*[left],
        count*[right],
        count*[fraction_low],
        count*[fraction_high],
        count*[quality_threshold],
        count*[max(verbose-1, 0)]
    )

    if (verbose):
        stderr.write('Calculating cleavage profile...\n')

    try:
        pool = Pool(workers, maxtasksperchild=500)
        # chunksize limited for memory
        interval_scores = pool.imap(
            _cleavage_profile_star,
            interval_list,
            chunksize=100
        )

        # output
        if (type(output_file) == str):   # check if output specified
            if (verbose):
                stderr.write(
                    f'Output file {output_file} specified. Opening...\n'
                )

            if output_file.endswith(".bw"): # BigWig
                with pbw.open(output_file, 'w') as bigwig:
                    bigwig.addHeader(header)
                    for interval_score in interval_scores:
                        contigs = interval_score['contig']
                        starts = interval_score['pos']
                        scores = interval_score['proportion']
                        stops = starts + 1

                        # skip empty intervals
                        if contigs.shape == (0,):
                            continue

                        try:
                            bigwig.addEntries(
                                chroms=contigs,
                                starts=starts,
                                ends=stops,
                                values=scores.astype(np.float64),
                            )
                        except RuntimeError:
                            stderr.write(
                                f'{contigs[0]}:{starts[0]}-{stops[-1]}\n'
                            )
                            stderr.write(
                                '/n invalid or out of order interval '
                                'encountered. Skipping to next.\n'
                            )
                            continue
            elif (output_file.endswith('.bed.gz')
                  or output_file.endswith('bedgraph.gz')
                  or output_file == "-"):
                  # XXX: writing to stdout is untested and may not work.
                with gzip.open(output_file, 'wt') as bedgraph:
                    for interval_score in interval_scores:
                        contigs = interval_score['contig']
                        starts = interval_score['pos']
                        scores = interval_score['proportion']
                        stops = starts + 1

                        lines = ''.join(f'{contig}\t{start}\t{stop}\t{score}\n'
                                 for contig, start, stop, score
                                 in zip(contigs, starts, stops, scores))

                        bedgraph.write(lines)

            else:   # unaccepted file type
                raise ValueError(
                    'output_file can only have suffix .bw'
                    )

        elif (output_file is not None):
            raise TypeError(
                f'output_file is unsupported type "{type(input_file)}". '
                'output_file should be a string specifying the path of the '
                'file to output scores to.'
                )
    finally:
        pool.close()

    if (verbose):
        end_time = time.time()
        stderr.write(
            f'cleavage profile took {end_time - start_time} s to complete\n'
        )
    
    return None