from __future__ import annotations
import gzip
import time
from multiprocessing.pool import Pool
from typing import Union

import pysam
import numpy as np
from numba import jit
from tqdm import tqdm

from finaletools.frag.wps import wps


def aggregate_wps(input_file: Union[pysam.AlignmentFile, str],
                  site_bed: str,
                  output_file: str=None,
                  window_size: int=120,
                  size_around_sites: int=5000,
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
        print(
            f"""
            Calculating aggregate WPS
            input_file: {input_file}
            site_bed: {site_bed}
            output_file: {output_file}
            window_size: {window_size}
            size_around_sites: {size_around_sites}
            quality_threshold: {quality_threshold}
            workers: {workers}
            verbose: {verbose}
            """
            )

    """
    with open(site_bed) as bed_file:
        contigs = []
        for line in bed_file:
            contig = line.split()[0].strip()
            if contig not in contigs:
                contigs.append(contig)

        num_contigs = len(contigs)


    if (verbose):
        print(f'Fragments for {num_contigs} contigs detected.')

    if (verbose >= 2):
        for contig in contigs:
            print(contig)

    input_tuples = zip([input_file] * num_contigs,
                       contigs,
                       [site_bed] * num_contigs,
                       [window_size] * num_contigs,
                       [size_around_sites] * num_contigs,
                       [fraction_low] * num_contigs,
                       [fraction_high] * num_contigs,
                       [quality_threshold] * num_contigs,
                       [verbose - 1 if verbose >= 1 else 0] * num_contigs)

    if (verbose):
        print('Calculating...')

    with Pool(workers) as pool:
        contig_scores = pool.starmap(_agg_wps_single_contig, input_tuples)

    if (verbose):
        print('Compiling scores')
    """
    # read tss contigs and coordinates from bed
    contigs = []
    starts = []
    stops = []
    with open(site_bed) as bed:
        for line in bed:
            contents = line.split()
            contig = contents[0].strip()
            start = int(contents[1])
            stop = int(contents[2])
            contigs.append(contig)
            starts.append(start)
            stops.append(stop)


    left_of_site = round(-size_around_sites / 2)
    right_of_site = round(size_around_sites / 2)

    assert right_of_site - left_of_site == size_around_sites

    count = len(contigs)

    tss_list = zip(
        count*[input_file],
        contigs,
        starts,
        stops,
        count*[None],
        count*[window_size],
        count*[fraction_low],
        count*[fraction_high],
        count*[quality_threshold])

    with Pool(workers) as pool:
        contig_scores = pool.starmap(wps, tss_list)

    scores = np.zeros((size_around_sites, 2))

    scores[:, 0] = np.arange(left_of_site, right_of_site)

    for contig_score in contig_scores:
        scores[:, 1] = scores[:, 1] + contig_score[:, 1]

    if (type(output_file) == str):   # check if output specified
        if (verbose):
            print(f'Output file {output_file} specified. Opening...')
        if output_file.endswith(".wig.gz"): # zipped wiggle
            with gzip.open(output_file, 'wt') as out:
                if (verbose):
                    print(f'File opened! Writing...')

                # declaration line
                out.write(
                    f'fixedStep\tchrom=.\tstart={left_of_site}\tstep={1}\tspan'
                    f'={window_size}\n'
                    )
                for score in (tqdm(scores[:, 1])
                              if verbose >= 2
                              else scores[:, 1]):
                    out.write(f'{score}\n')

        elif output_file.endswith(".wig"):  # wiggle
            with open(output_file, 'wt') as out:
                if (verbose):
                    print(f'File opened! Writing...')
                # declaration line
                out.write(
                    f'fixedStep\tchrom=.\tstart={left_of_site}\tstep={1}\tspan'
                    f'={window_size}\n'
                    )
                for score in (tqdm(scores[:, 1])
                              if verbose >= 2
                              else scores[:, 1]):
                    out.write(f'{score}\n')

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
        print(f'aggregate_wps took {end_time - start_time} s to complete',
              flush=True)

    return scores