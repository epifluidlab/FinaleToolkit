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


def delfi(input_bam: Union[str, pysam.AlignmentFile],
          genome_file: str,
          blacklist_file: str,
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
    input_bam: str or AlignmentFile
        Path string or Alignment File pointing a bam file containing PE
        fragment reads.
    genome_file: str
        Path string to .genome file.
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

    with open(genome_file) as genome:
        # read genome file into a list of tuples
        contigs = [(
            line.split()[0],
            int(line.split()[1])
            ) for line in genome.readlines()]

    # generate DELFI windows
    windows = []
    for contig, size in contigs:
        for coordinate in range(0, size, window_size):
            # (contig, start, stop)
            windows.append((contig, coordinate, coordinate + window_size))

    with Pool(workers) as pool:
        pass

    print(contigs)



    if (verbose):
        end_time = time.time()
        print(f'aggregate_wps took {end_time - start_time} s to complete')
    return None
