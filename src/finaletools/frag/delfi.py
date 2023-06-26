from __future__ import annotations
import argparse
import gzip
import time
import os
import tempfile as tf
from multiprocessing.pool import Pool
from typing import Union, TextIO, BinaryIO

import pysam
import numpy as np
from numba import jit
from tqdm import tqdm
from finaletools.utils import (
    frag_bam_to_bed,
    frag_array,
    not_read1_or_low_quality
)


def _delfi_single_window(
        input_file: str,
        contig: str,
        window_start: int,
        window_stop: int,
        blacklist_file: str=None,
        verbose: Union[int,bool]=False) -> int:
    """
    Calculates DELFI for one window.
    """

    if (blacklist_file is not None):
        blacklist_regions = []

        # convert blacklist to a list of tuples
        # TODO: accept other file types
        with open(blacklist_file) as blacklist:
            for line in blacklist:
                contents = line.split()
                blacklist_regions.append(contents[0],
                                         int(contents[1]),
                                         int(contents[2])
                )

        


    return None


def delfi(input_file: str,  # TODO: allow AlignmentFile to be used
          genome_file: str,
          blacklist_file: str=None,
          window_size: int=5000000,
          subsample_coverage: int=2,
          quality_threshold: int=30,
          workers: int=1,
          preprocessing: bool=True,
          verbose: Union[int, bool]=False):
    """
    A function that replicates the methodology of Christiano et al
    (2019).

    Parameters
    ----------
    input_file: str
        Path string pointing to a bam file containing PE
        fragment reads.
    genome_file: str
        Path string to .genome file. Should contain only autosomal chromosomes
    blacklist_file: str
        Path string to bed file containing genome blacklist.
    window_size: int
        Size of non-overlapping windows to cover genome. Default is
        5 megabases.
    subsample_coverage: int, optional
        The depth at which to subsample the input_bam. Default is 2.
    workers: int, optional
        Number of worker processes to use. Default is 1.
    preprocessing: bool, optional
        Christiano et al (2019)
    verbose: int or bool, optional
        Determines how many print statements and loading bars appear in
        stdout. Default is False.

    """

    # TODO: subsample bam to specified coverage. Jan28.hg19.mdups.bam
    # already has 1-2x coverage.

    if (verbose):
        start_time = time.time()

    contigs = []
    with open(genome_file) as genome:
        for line in genome:
            contents = line.split()
            contigs.append((contents[0],  int(contents[1])))


    # generate DELFI windows
    windows = []
    for contig, size in contigs:
        for coordinate in range(0, size, window_size):
            # (contig, start, stop)
            windows.append((input_file,
                            contig,
                            coordinate,
                            coordinate + window_size,
                            blacklist_file))

    with Pool(workers) as pool:
        pass

    print(contigs)



    if (verbose):
        end_time = time.time()
        print(f'aggregate_wps took {end_time - start_time} s to complete')
    return None
