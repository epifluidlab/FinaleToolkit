from __future__ import annotations
import gzip
import time
from multiprocessing.pool import Pool
from typing import Union
from sys import stderr, stdout

import pysam
import numpy as np
from numba import jit
from tqdm import tqdm
from memory_profiler import profile

from finaletools.frag.wps import wps


def multi_wps(input_file: Union[pysam.AlignmentFile, str],
                  site_bed: str,
                  output_file: str=None,
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
        BAM or SAM file containing paired-end fragment reads or its
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

    # read tss contigs and coordinates from bed
    if (verbose):
        stderr.write('Reading intervals from bed\n')

    contigs = []
    starts = []
    stops = []
    strands = []
    with open(site_bed) as bed:
        for line in bed:
            contents = line.split()
            contig = contents[0].strip()
            start = int(contents[1])
            stop = int(contents[2])
            strand = contents[5].strip()
            contigs.append(contig)
            starts.append(start)
            stops.append(stop)
            strands.append(strand)


    left_of_site = round(-interval_size / 2)
    right_of_site = round(interval_size / 2)

    assert right_of_site - left_of_site == interval_size

    count = len(contigs)

    if (verbose):
        stderr.write('Zipping inputs\n')

    tss_list = zip(
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

    with Pool(workers, maxtasksperchild=500) as pool:
        contig_scores = pool.starmap(wps, tss_list, chunksize=10000)

    scores = np.zeros((interval_size, 2))

    scores[:, 0] = np.arange(left_of_site, right_of_site)

    if (verbose):
        stderr.write('Flipping scores for reverse strand\n')

    # aggregate scores by strand
    for i in range(count):
        if strands[i] == '.':
            continue
        elif strands[i] == '-':
            contig_score = np.flip(contig_scores[i], 1)
            scores[:, 1] = scores[:, 1] + contig_score[:, 1]
        elif strands[i] == '+':
            contig_score = contig_scores[i]
            scores[:, 1] = scores[:, 1] + contig_score[:, 1]
        else:
            stderr.write('Invalid strand found. Interval skipped.')


    if (type(output_file) == str):   # check if output specified
        if (verbose):
            stderr.write(f'Output file {output_file} specified. Opening...\n')
        if output_file.endswith(".wig.gz"): # zipped wiggle
            with gzip.open(output_file, 'wt') as out:
                if (verbose):
                    stderr.write(f'File opened! Writing...\n')

                # declaration line
                out.write(
                    f'fixedStep\tchrom=.\tstart={left_of_site}\tstep={1}\tspan'
                    f'={window_size}\n'
                    )
                for score in (tqdm(scores[:, 1]+1)
                              if verbose >= 2
                              else scores[:, 1]+1):
                    out.write(f'{score}\n')

        elif output_file.endswith(".wig"):  # wiggle
            with open(output_file, 'wt') as out:
                if (verbose):
                    stderr.write(f'File opened! Writing...\n')
                # declaration line
                out.write(
                    f'fixedStep\tchrom=.\tstart={left_of_site}\tstep={1}\tspan'
                    f'={interval_size}\n'
                    )
                for score in (tqdm(scores[:, 1]+1)
                              if verbose >= 2
                              else scores[:, 1]+1):
                    out.write(f'{score}\n')

        elif output_file == '-':  # stdout
            if (verbose):
                stderr.write(f'File opened! Writing...\n')
            # declaration line
            stdout.write(
                f'fixedStep\tchrom=.\tstart={left_of_site}\tstep={1}\tspan'
                f'={window_size}\n'
                )
            for score in (tqdm(scores[:, 1])
                            if verbose >= 2
                            else scores[:, 1]):
                stdout.write(f'{score}\n')

        else:   # unaccepted file type
            raise ValueError(
                'output_file can only have suffixes .wig or .wig.gz.'
                )

    elif (output_file is not None):
        raise TypeError(
            f'output_file is unsupported type "{type(input_file)}". '
            'output_file should be a string specifying the path of the file '
            'to output scores to.'
            )

    if (verbose):
        end_time = time.time()
        stderr.write(
            f'aggregate_wps took {end_time - start_time} s to complete\n'
        )

    return scores