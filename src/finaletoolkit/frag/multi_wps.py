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

from finaletoolkit.frag.wps import wps


def _wps_star(args):
    """Helper function to be used with imap"""
    return wps(*args)


def multi_wps(
        input_file: Union[pysam.AlignmentFile, str],
        site_bed: str,
        output_file: Union[str, None]=None,
        window_size: int=120,
        interval_size: int=5000,
        fraction_low: int=120,
        fraction_high: int=180,
        quality_threshold: int=30,
        workers: int=1,
        verbose: Union[bool, int]=0
        ) -> np.ndarray:
    """
    Function that aggregates WPS over sites in BED file according to the
    method described by Snyder et al (2016).

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, SAM, or tabix file containing paired-end fragment reads or its
        path. `AlignmentFile` must be opened in read mode.
    site_bed: str
        Bed file containing intervals to perform WPS on.
    output_file : string, optional
    window_size : int, optional
        Size of window to calculate WPS. Default is k = 120, equivalent
        to L-WPS.
    interval_size : int, optional
        Size of each interval specified in the bed file. Should be the
        same for every interbal. Default is 5000.
    fraction_low : int, optional
        Specifies lowest fragment length included in calculation.
        Default is 120, equivalent to long fraction.
    fraction_high : int, optional
        Specifies highest fragment length included in calculation.
        Default is 120, equivalent to long fraction.
    quality_threshold : int, optional
    workers : int, optional
    verbose : bool, optional

    Returns
    -------
    scores : numpy.ndarray
        np array of shape (n, 2) where column 1 is the coordinate and
        column 2 is the score and n is the number of coordinates in
        region [start,stop)
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
        if site_bed != '-':
            bed.close()

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
        count*[fraction_low],
        count*[fraction_high],
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

    return scores
