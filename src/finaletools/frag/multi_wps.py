from __future__ import annotations
import gzip
import time
from multiprocessing.pool import Pool
from typing import Union, Iterator, Generator
from sys import stderr, stdout, stdin

import pysam
import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm
import pyBigWig as pbw

from finaletools.frag.wps import wps
from finaletools.utils.utils import _parse_chrom_sizes, _reduce_overlaps_in_file, _convert_to_list, _merge_all_intervals

def _wps_star(args):
    return wps(*args)


def multi_wps(
        input_file: Union[pysam.AlignmentFile, str],
        interval_file: str,
        output_file: Union[str, None],
        window_size: int,
        interval_size: int,
        min_length: int,
        max_length: int,
        chrom_sizes: str,
        quality_threshold: int,
        intersect_policy: str,
        workers: int,
        verbose: bool
        ) -> np.ndarray:
   
    if (verbose):
        start_time = time.time()
        stderr.write(
            f"""
            Calculating aggregate WPS
            input_file: {input_file}
            interval_file: {interval_file}
            output_file: {output_file}
            window_size: {window_size}
            min_length: {min_length}
            max_length: {max_length}
            chrom_sizes: {chrom_sizes}
            interval_size: {interval_size}
            intersect_policy: {intersect_policy}
            quality_threshold: {quality_threshold}
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

    if (verbose):
        stderr.write('Reading intervals from bed\n')

    reduced_intervals = _reduce_overlaps_in_file(interval_file) # WARNING THIS NEEDS TO BE CHANGED BECAUSE IT DOES NOT LEAVE ADJACENT INTERVALS INTACT
    converted_intervals = _convert_to_list(reduced_intervals)
    all_intervals = _merge_all_intervals(converted_intervals)
    contigs, starts, stops = zip(*all_intervals)
    count = len(contigs)

    left_of_site = round(-interval_size / 2)
    right_of_site = round(interval_size / 2)

    assert right_of_site - left_of_site == interval_size

    count = len(contigs)

    if (verbose):
        stderr.write('Zipping inputs\n')

    interval_list = zip(
        count*[input_file],
        contigs,
        starts,
        stops,
        count*[None],
        count*[window_size],
        count*[min_length],
        count*[max_length],
        count*[intersect_policy],
        count*[quality_threshold],
        count*[verbose]
    )

    if (verbose):
        stderr.write('Calculating WPS...\n')

    try:
        pool = Pool(processes=workers)
        interval_scores = pool.imap(_wps_star, interval_list, chunksize=max(len(contigs)//workers, 1))

        if (type(output_file) == str):
            if (verbose):
                stderr.write(
                    f'Output file {output_file} specified. Opening...\n'
                )

            if output_file.endswith(".bw"): # BigWig
                with pbw.open(output_file, 'w') as bigwig:
                    bigwig.addHeader(header)
                    for interval_score in tqdm(interval_scores, desc="Adding intervals to bigWig file.", disable=not verbose):
                        contigs = interval_score['contig']
                        starts = interval_score['start']
                        scores = interval_score['wps']
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
            else:
                raise ValueError(
                    'Output file can only have suffix .bw'
                    )

        elif (output_file is not None):
            raise TypeError(
                f'The output file is unsupported type "{type(output_file)}". '
                'Output file should be a string specifying the path of the '
                'file to output scores to.'
                )
    finally:
        pool.close()

    if (verbose):
        end_time = time.time()
        stderr.write(
            f'Calculating WPS took {end_time - start_time} s to complete\n'
        )

    return scores
