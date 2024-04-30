"""Cleavage Profiler

This module is used to find the cleavage profile as described in Zhou
et al 2022 (https://doi.org/10.1073/pnas.2209852119). Cleavage profile
describes the proportion of fragment ends at a site over the depth at
the site (cleavage proportion) calculated over a 5+/- nt window around a
CpG site.
"""

from __future__ import annotations
from typing import Union
from sys import stderr, stdout, stdin
from multiprocessing import Pool
import gzip
import time

import numpy as np
import pyBigWig as pbw
from pybedtools import BedTool
import pysam

from finaletoolkit.utils.utils import frag_array, overlaps


def cleavage_profile(
    input_file: str,
    contig: str,
    start: int,
    stop: int,
    fraction_low: int=1,
    fraction_high: int=10000000,
    quality_threshold: int=30,
    verbose: Union[bool, int]=0
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
            fraction_low: {fraction_low}
            fraction_high: {fraction_high}
            quality_threshold: {quality_threshold}
            verbose: {verbose}
            """
        )

    frags = frag_array(
        input_file=input_file,
        contig=contig,
        quality_threshold=quality_threshold,
        start=start,
        stop=stop,
        fraction_low=fraction_low,
        fraction_high=fraction_high
    )

    positions = np.arange(start, stop)

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
    proportions = ends/depth*100

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


def _cli_cleavage_profile(
    input_file: str,
    interval_file: str,
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
    
      # get header from input_file
    if (input_file.endswith('.sam')
        or input_file.endswith('.bam')
        or input_file.endswith('.cram')):
        with pysam.AlignmentFile(input_file, 'r') as bam:
            references = bam.references
            lengths = bam.lengths
            header = list(zip(references, lengths))
    # TODO: get a header when reading tabix
    elif (input_file.endswith('.bed')
          or input_file.endswith('.bed.gz')
          or input_file.endswith('.frag')
          or input_file.endswith('.frag.gz')
    ):
        with pysam.TabixFile(input_file, 'r') as tbx:
            raise NotImplementedError('tabix files not yet supported!')
    else:
        raise ValueError("Not a supported file type.")

    if (verbose > 1):
        stderr.write(f'header is {header}\n')
    
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
            contents = line.split()
            contig = contents[0].strip()
            start = int(contents[1])
            stop = int(contents[2])

            # cut off part of previous interval if overlap
            if contig == prev_contig and start < prev_stop:
                prev_stop = start

            if prev_contig is not None:
                contigs.append(prev_contig)
                starts.append(prev_start)
                stops.append(prev_stop)

            prev_contig = contig
            prev_start = start
            prev_stop = stop
        # appending last interval
        contigs.append(prev_contig)
        starts.append(prev_start)
        stops.append(prev_stop)
    finally:
        if interval_file != '-':
            bed.close()
    
    count = len(contigs)

    if (verbose):
        stderr.write('Zipping inputs\n')

    interval_list = zip(
        count*[input_file],
        contigs,
        starts,
        stops,
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
                        except RuntimeError as e:
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