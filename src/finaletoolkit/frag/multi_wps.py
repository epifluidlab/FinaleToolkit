from __future__ import annotations
import gzip
import time
from multiprocessing.pool import Pool
from typing import Union
from sys import stderr, stdin
import warnings

import pysam
import numpy as np
import pyBigWig as pbw

from finaletoolkit.frag.wps import wps
from finaletoolkit.utils.utils import chrom_sizes_to_list


def _wps_star(args):
    """Helper function to be used with imap"""
    return wps(*args)


def multi_wps(
        input_file: Union[pysam.AlignmentFile, str],
        site_bed: str,
        chrom_sizes: str=None,
        output_file: Union[str, None]=None,
        window_size: int=120,
        interval_size: int=5000,
        min_length: int=120,
        max_length: int=180,
        quality_threshold: int=30,
        workers: int=1,
        verbose: Union[bool, int]=0,
        fraction_low: int=None,
        fraction_high: int=None,
        ):
    """
    Function that aggregates WPS over sites in BED file according to the
    method described by Snyder et al (2016).

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, SAM, or tabix file containing paired-end fragment reads or its
        path. `AlignmentFile` must be opened in read mode.
    site_bed: str
        BED file containing intervals to perform WPS on. The intervals
        in this BED file should be sorted, first by `contig` then
        `start`.
    chrom_sizes: str or pathlike, optional
        Tab separated file containing names and sizes of chromosomes in
        `input_file`. Required if `input_file` is tabix-indexed.
    output_file : string, optional
    window_size : int, optional
        Size of window to calculate WPS. Default is k = 120, equivalent
        to L-WPS.
    interval_size : int, optional
        Size of each interval specified in the bed file. Should be the
        same for every interval. Default is 5000.
    min_length : int, optional
        Specifies lowest fragment length included in calculation.
        Default is 120, equivalent to long fraction.
    max_length : int, optional
        Specifies highest fragment length included in calculation.
        Default is 120, equivalent to long fraction.
    quality_threshold : int, optional
    workers : int, optional
    verbose : bool, optional
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
            Calculating aggregate WPS
            input_file: {input_file}
            site_bed: {site_bed}
            output_file: {output_file}
            window_size: {window_size}
            interval_size: {interval_size}
            quality_threshold: {quality_threshold}
            workers: {workers}
            verbose: {verbose}

            """
        )

    if (input_file == '-' and site_bed == '-'):
        raise ValueError('input_file and site_bed cannot both read from stdin')
    
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

    # get chrom sizes from input_file or chrom_sizes
    if (input_file.endswith('.sam')
        or input_file.endswith('.bam')
        or input_file.endswith('.cram')):
        with pysam.AlignmentFile(input_file, 'r') as bam:
            references = bam.references
            lengths = bam.lengths
            header = list(zip(references, lengths))
    elif (isinstance(input_file, pysam.AlignmentFile)):
        pass
    else:
        if chrom_sizes is None:
            raise ValueError(
                'chrom_sizes must be specified for BED/Fragment files'
            )
        header = chrom_sizes_to_list(chrom_sizes)

    if (verbose > 1):
        stderr.write(f'header is {header}\n')

    # read tss contigs and coordinates from bed
    if (verbose):
        stderr.write('Reading intervals from bed\n')

    contigs = []
    starts = []
    stops = []
    try:
        if site_bed == '-':
            bed = stdin
        else:
            bed = open(site_bed)

        left_of_site = round(-interval_size / 2)
        right_of_site = round(interval_size / 2)

        assert right_of_site - left_of_site == interval_size

        # for overlap checking
        prev_contig = None
        prev_start = 0
        prev_stop = 0
        for line in bed:
            contents = line.split()
            contig = contents[0].strip()
            midpoint = (int(contents[1]) + int(contents[2])) // 2

            start = max(0, midpoint + int(left_of_site))
            stop = midpoint + int(right_of_site)

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
        if site_bed != '-':
            bed.close()  

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
        count*[quality_threshold],
        count*[verbose-2 if verbose>2 else 0]
    )

    if (verbose):
        stderr.write('Calculating wps...\n')

    try:
        pool = Pool(workers, maxtasksperchild=500)
        # chunksize limited for memory
        interval_scores = pool.imap(
            _wps_star,
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
                        starts = interval_score['start']
                        scores = interval_score['wps']

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
                  or output_file.endswith('bedGraph.gz')):
                with gzip.open(output_file, 'wt') as bedgraph:
                    for interval_score in interval_scores:
                        contigs = interval_score['contig']
                        starts = interval_score['start']
                        scores = interval_score['wps']
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
            f'multi_wps took {end_time - start_time} s to complete\n'
        )
    return output_file
