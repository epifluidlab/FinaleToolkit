# TODO: still wip, need to fully implement

from __future__ import annotations
import time
from multiprocessing.pool import Pool
from typing import Union
from tempfile import TemporaryDirectory

import pysam
import py2bit
import numpy as np
from finaletools.utils import not_read1_or_low_quality


def _delfi_single_window(
        input_file: str,
        reference_file: str,
        contig: str,
        window_start: int,
        window_stop: int,
        blacklist_file: str=None,
        quality_threshold: int=30,
        verbose: Union[int,bool]=False) -> tuple:
    """
    Calculates DELFI for one window.
    """

    blacklist_regions = []

    if (blacklist_file is not None):


        # convert blacklist to a list of tuples
        # TODO: accept other file types
        with open(blacklist_file) as blacklist:
            for line in blacklist:
                region_contig, region_start, region_stop, *_ = line.split()
                region_start = int(region_start)
                region_stop = int(region_stop)
                if (contig == region_contig
                    and window_start <= region_start
                    and window_stop >= region_stop ):
                    blacklist_regions.append((region_start,region_stop))

    short_lengths = []
    long_lengths = []
    gc_tally = 0  # cumulative sum of gc bases
    base_tally = 0

    num_frags = 0


    with pysam.AlignmentFile(input_file) as sam_file:
        # Iterating on each read in file in specified contig/chromosome
        for read1 in (sam_file.fetch(contig=contig,
                                        start=window_start,
                                        stop=window_stop)):

            # Only select forward strand and filter out non-paired-end
            # reads and low-quality reads
            if (not_read1_or_low_quality(read1, quality_threshold)):
                pass
            else:
                frag_start = read1.reference_start
                frag_length = read1.template_length
                frag_end = frag_start + frag_length

                # check if in blacklist
                blacklisted = False
                for region in blacklist_regions:
                    if (
                        (frag_start >= region[0] and frag_start < region[1])
                        and (frag_end >= region[0] and frag_end < region[1])
                    ):
                        blacklisted = True
                        break


                if (not blacklisted
                    and frag_length >= 100
                    and frag_length <= 220
                ):
                    # append length of fragment to list
                    if (frag_length >= 151):
                        long_lengths.append(abs(frag_length))
                    else:
                        short_lengths.append(abs(frag_length))

                    num_frags += 1


    with py2bit.open(reference_file) as ref_seq:
        ref_bases = ref_seq.sequence(contig, window_start, window_stop).upper()

    num_bases = window_stop - window_start
    num_gc = sum([(base == 'G' or base == 'C') for base in ref_bases])

    # NaN if no fragments in window.
    gc_content = num_gc / num_bases if num_frags > 0 else np.NaN

    coverage_short = len(short_lengths)
    coverage_long = len(long_lengths)

    # if (len(short_lengths) != 0 or len(long_lengths) != 0):
    if verbose:
        print(f'{contig}:{window_start}-{window_stop} short: '
              f'{coverage_short} long: {coverage_long}, gc_content: '
              f'{gc_content*100}%')

    return (contig,
            window_start,
            window_stop,
            coverage_short,
            coverage_long,
            gc_content,
            num_frags)


def delfi(input_file: str,  # TODO: allow AlignmentFile to be used
          autosomes: str,
          reference_file: str,
          blacklist_file: str=None,
          output_file: str=None,
          window_size: int=100000,
          subsample_coverage: float=2,
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
    autosomes: str
        Path string to a .genome file containing only autosomal chromosomes
    reference_file: str
        Path string to .2bit file.
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

    # TODO: add support to fasta for reference_file

    # TODO: add option to supply a reference genome name instead of
    # genome file

    # TODO: subsample bam to specified coverage. Jan28.hg19.mdups.bam
    # already has 1-2x coverage.

    if (verbose):
        start_time = time.time()
        print(f"""
        input_file: {input_file}
        autosomes: {autosomes}
        reference_file: {reference_file}
        blacklist_file: {blacklist_file}
        output_file: {output_file}
        window_size: {window_size}
        subsample_coverage: {subsample_coverage}
        quality_threshold: {quality_threshold}
        workers: {workers}
        preprocessing: {preprocessing}
        verbose: {verbose}
        """)

    if verbose:
        print(f'Reading genome file...')

    contigs = []
    with open(autosomes) as genome:
        for line in genome:
            contents = line.split('\t')
            # account for empty lines
            if len(contents) > 1:
                contigs.append((contents[0],  int(contents[1])))

    if verbose:
        print(f'Generating windows')

    # generate DELFI windows
    window_args = []
    for contig, size in contigs:
        for coordinate in range(0, size, window_size):
            # (contig, start, stop)
            window_args.append((input_file,
                                reference_file,
                                contig,
                                coordinate,
                                coordinate + window_size,
                                blacklist_file,
                                quality_threshold,
                                verbose - 1 if verbose > 1 else 0))

    if (verbose):
        print(f'{len(window_args)} windows created.')
        print('Calculating fragment lengths...')

    with Pool(workers) as pool:
        windows = pool.starmap(_delfi_single_window, window_args)

    with open(output_file, 'w') as out:
        out.write('contig\tstart\tstop\tshort\tlong\tgc%\n')
        for window in windows:
            out.write(f'{window[0]}\t{window[1]}\t{window[2]}\t{window[3]}\t{window[4]}\t{window[5]}\n')


    num_frags = sum(window[3] for window in windows)

    if (verbose):
        end_time = time.time()
        print(f'{num_frags} fragments included.')
        print(f'delfi took {end_time - start_time} s to complete')
    return None
