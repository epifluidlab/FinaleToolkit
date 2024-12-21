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
from os import PathLike
from sys import stderr, stdin
from multiprocessing import Pool
import gzip
import time
import warnings

import numpy as np
import pysam
import pyBigWig as pbw

from finaletoolkit.utils.utils import (
    frag_array, chrom_sizes_to_list, _reduce_overlaps_in_file,
    _convert_to_list, _merge_all_intervals, chrom_sizes_to_dict
    )
from finaletoolkit.utils.typing import FragFile


def cleavage_profile(
    input_file: FragFile,
    chrom_size: int,
    contig: str,
    start: int,
    stop: int,
    left: int=0,
    right: int=0,
    min_length: int|None=None,
    max_length: int|None=None,
    quality_threshold: int=30,
    verbose: Union[bool, int]=0,
    fraction_low: int|None=None,
    fraction_high: int|None=None,
) -> np.ndarray:
    """
    Cleavage profile calculated over a single interval.

    Parameters
    ---------
    input_file: str
        BAM, CRAM, or FRAG file with fragment information.
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
    min_length: int
        Minimum fragment size to include
    max_length: int
        Maximum fragment size to include
    quality_threshold: int
        Minimum MAPQ
    verbose: bool or in
    fraction_low : int, optional
        Deprecated alias for min_length
    fraction_high : int, optional
        Deprecated alias for max_length

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

    # Pass aliases and check for conflicts
    if fraction_low is not None and min_length is None:
        min_length = fraction_low
        warnings.warn("fraction_low is deprecated. Use min_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
    elif fraction_low is not None and min_length is not None:
        warnings.warn("fraction_low is deprecated. Use min_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
        raise ValueError(
            'fraction_low and min_length cannot both be specified')

    if fraction_high is not None and max_length is None:
        max_length = fraction_high
        warnings.warn("fraction_high is deprecated. Use max_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
    elif fraction_high is not None and max_length is not None:
        warnings.warn("fraction_high is deprecated. Use max_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
        raise ValueError(
            'fraction_high and max_length cannot both be specified')

    adj_start = max(start-left, 0)
    adj_stop = min(stop+right, chrom_size)

    frags = frag_array(
        input_file=input_file,
        contig=contig,
        quality_threshold=quality_threshold,
        start=adj_start,
        stop=adj_stop,
        min_length=min_length,
        max_length=max_length,
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


def multi_cleavage_profile(
    input_file: FragFile,
    interval_file: Union[str, Path],
    chrom_sizes: Union[str, Path, None] = None,
    left: int=0,
    right: int=0,
    min_length: int|None=None,
    max_length: int|None=None,
    quality_threshold: int=30,
    output_file: str='-',
    workers: int=1,
    verbose: Union[bool, int]=0,
    fraction_low: int|None=None,
    fraction_high: int|None=None,
):
    """
    Multithreaded cleavage profile implementation over intervals in a
    BED file.
    Parameters
    ---------
    input_file: str or pathlike
        BAM, CRAM, or FRAG file with fragment information.
    interval_file: str or pathlike
        Sorted BED file containing intervals.
    chrom_sizes: str or pathlike
        Tab-delimited file with name and lengths of each contig. Required
        if input_file is a tabix-indexed frag file.
    left: int
        Amount to subtract from start coordinate. Useful if only given
        coordinates of CpG.
    right: int
        Amount to add to stop coordinate.
    min_length: int
        Minimum fragment size to include
    max_length: int
        Maximum fragment size to include
    quality_threshold: int
        Minimum MAPQ
    output_file: str or Pathlike
        Bigwig or bedgraph file to write results to.
    workers: int, default = 1
        Number of processes to spawn
    verbose: bool or int
    fraction_low : int, optional
        Deprecated alias for min_length
    fraction_high : int, optional
        Deprecated alias for max_length
        
    Returns
    -------
    output_file: str
        location results are stored.
    """

    if (verbose):
        start_time = time.time()
        stderr.write(
            f"""
            Calculating cleavage profile
            input_file: {input_file}
            interval_file: {interval_file}
            chrom_sizes: {chrom_sizes}
            left: {left}
            right: {right}
            min_length: {min_length}
            max_length: {max_length}
            quality_threshold: {quality_threshold}
            output_file: {output_file}
            workers: {workers}
            verbose: {verbose}
            """
        )
    # Pass aliases and check for conflicts
    if fraction_low is not None and min_length is None:
        min_length = fraction_low
        warnings.warn("fraction_low is deprecated. Use min_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
    elif fraction_low is not None and min_length is not None:
        warnings.warn("fraction_low is deprecated. Use min_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
        raise ValueError(
            'fraction_low and min_length cannot both be specified')

    if fraction_high is not None and max_length is None:
        max_length = fraction_high
        warnings.warn("fraction_high is deprecated. Use max_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
    elif fraction_high is not None and max_length is not None:
        warnings.warn("fraction_high is deprecated. Use max_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
        raise ValueError(
            'fraction_high and max_length cannot both be specified')
    
    if (input_file == '-' and interval_file == '-'):
        raise ValueError('input_file and site_bed cannot both read from stdin')
    
    if chrom_sizes is None:
        raise ValueError(
            '--chrom_sizes must be specified'
        )
    # get chroms
    header = chrom_sizes_to_list(chrom_sizes)
    chrom_dict = chrom_sizes_to_dict(chrom_sizes)

    if (verbose > 1):
        stderr.write(f'chrom sizes {header}\n')
    
    # reading intervals from bed and removing overlaps
    # NOTE: assumes that bed file is sorted.
    contigs = []
    starts = []
    stops = []
    try:
        if interval_file == '-':
            bed = stdin
        else:
            bed = open(interval_file)

        # for overlap checking
        prev_contig = None
        prev_start = 0
        prev_stop = 0
        for line in bed:
            # parse file
            contents = line.split()
            contig = contents[0].strip()
            start, stop = int(contents[1]), int(contents[2])
            start = max(0, start - left)
            stop = min(stop + right, chrom_dict[contig])

            # if overlapping:
            if (prev_contig == contig and start < prev_stop):
                prev_stop = max(prev_stop, stop)
            else:
                contigs.append(prev_contig)
                starts.append(prev_start)
                stops.append(prev_stop)
                prev_contig, prev_start, prev_stop = contig, start, stop
        contigs.append(prev_contig)
        starts.append(prev_start)
        stops.append(prev_stop)
    finally:
        if interval_file != '-':
            bed.close()  
            
    contigs = contigs[1:]
    starts = starts[1:]
    stops = stops[1:]

    # reading chrom.sizes file
    # get chrom sizes from input_file or chrom_sizes
    
    if (isinstance(input_file, pysam.AlignmentFile)):
        references = input_file.references
        lengths = input_file.lengths
        header = list(zip(references, lengths))
    elif (
            (isinstance(input_file, str)
            or isinstance(input_file, PathLike))
        and
            (str(input_file).endswith('.sam')
            or str(input_file).endswith('.bam')
            or str(input_file).endswith('.cram'))
    ):
        with pysam.AlignmentFile(input_file, 'r') as bam:
            references = bam.references
            lengths = bam.lengths
            header = list(zip(references, lengths))
    elif chrom_sizes is not None:
        header = chrom_sizes_to_list(chrom_sizes)
    else:
        raise ValueError(
            'chrom_sizes must be specified for Tabix-indexed files'
        )
        
    size_dict = dict(header)

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
        count*[0],  # left and right precomputed to avoid overlaps
        count*[0],
        count*[min_length],
        count*[max_length],
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
                    last = "None"
                    for interval_score in interval_scores:
                        
                        contigs = interval_score['contig']
                        starts = interval_score['pos']
                        scores = interval_score['proportion']

                        # skip empty intervals
                        if contigs.shape == (0,):
                            continue

                        try:

                            bigwig.addEntries(
                                contigs[0],
                                starts[0],
                                values=scores.astype(np.float64),
                                step=1,
                                span=1
                            )
                        except RuntimeError as e:
                            stderr.write(
                                f'{contigs[0]}:{starts[0]}-{starts[-1]+1}\n'
                            )
                            stderr.write(
                                'invalid or out of order interval '
                                'encountered. Skipping to next.\n'
                            )
                            stderr.write(
                                f"captured error:\n{e}\n")
                            stderr.write(
                                f"current output:\n{interval_score}\n")
                            stderr.write(
                                f"last output:\n{last}\n")
                            continue
                        last = interval_score
            elif (output_file.endswith('.bed.gz')
                  or output_file.endswith('bedgraph.gz')
                  or output_file == "-"):
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
                    'output_file can only have suffix .bw, .bedgraph.gz, or .bed.gz.'
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
    return output_file