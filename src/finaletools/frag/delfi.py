# TODO: still wip, need to fully implement

from __future__ import annotations
import time
from multiprocessing.pool import Pool
from typing import Union, TextIO
from tempfile import TemporaryDirectory
from sys import stderr, stdout
import gzip

import pysam
import py2bit
import numpy as np

from finaletools.utils.utils import _not_read1_or_low_quality


def _delfi_single_window(
        input_file: str,
        reference_file: str,
        contig: str,
        window_start: int,
        window_stop: int,
        blacklist_file: str=None,
        tcmeres: list=None,
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
    frag_pos = []

    num_frags = 0

    try:
        # read from tabix or bam/bam
        if input_file.endswith('.bam') or input_file.endswith('.sam'):
            file = pysam.AlignmentFile(input_file)
            is_sam = True
        elif (
            input_file.endswith('.bed')
            or input_file.endswith('.bed.gz')
            or input_file.endswith('.frag.gz')
            or input_file.endswith('.frag')
        ):
            file = pysam.TabixFile(input_file, parser=pysam.asBed())
            is_sam = False
        else:
            raise ValueError(
                'Unsupported type. Only BAM, SAM, and tabix indexed files'
                'accepted.'
            )
        # Iterating on each read in file in specified contig/chromosome
        for read1 in (file.fetch(contig, window_start, window_stop)):

            # Only select forward strand and filter out non-paired-end
            # reads and low-quality reads.
            # relies on short circuit evaluation to avoid calling sam
            # method on tabix.
            if is_sam and _not_read1_or_low_quality(read1, quality_threshold):
                pass
            else:
                # TODO: fix tabix reading
                if is_sam:
                    frag_start = read1.reference_start
                    frag_length = read1.template_length
                    frag_stop = frag_start + frag_length
                else:
                    frag_start = read1.start
                    frag_stop = read1.end
                    frag_length = frag_stop - frag_start

                # check if in blacklist
                blacklisted = False
                for region in blacklist_regions:
                    if (
                        (frag_start >= region[0] and frag_start < region[1])
                        and (frag_stop >= region[0] and frag_stop < region[1])
                    ):
                        blacklisted = True
                        break

                # check if in centromere or telomere
                in_tcmere = False
                for tc_contig, tc_start, tc_stop, _ in tcmeres:
                    if (tc_contig == region_contig
                        and (frag_start >= tc_start and frag_start < tc_stop)
                        and (frag_stop >= tc_start and frag_stop < tc_stop)
                    ):
                        in_tcmere = True
                        break

                if (not blacklisted
                    and not in_tcmere
                    and frag_length >= 100
                    and frag_length <= 220
                ):
                    # append length of fragment to list
                    if (frag_length >= 151):
                        long_lengths.append(abs(frag_length))
                    else:
                        short_lengths.append(abs(frag_length))

                    frag_pos.append((frag_start, frag_stop))

                    num_frags += 1
    finally:
        file.close()

    num_gc = 0  # cumulative sum of gc bases

    with py2bit.open(reference_file) as ref_seq:
        ref_bases = ref_seq.sequence(contig, window_start, window_stop).upper()

    num_gc = sum((base == 'G' or base == 'C') for base in ref_bases)

    # window_length
    window_coverage = window_stop - window_start

    # NaN if no fragments in window.
    gc_content = num_gc / window_coverage if num_frags > 0 else np.NaN

    coverage_short = len(short_lengths)
    coverage_long = len(long_lengths)

    # if (len(short_lengths) != 0 or len(long_lengths) != 0):
    if verbose:
        stderr.write(
            f'{contig}:{window_start}-{window_stop} short: '
            f'{coverage_short} long: {coverage_long}, gc_content: '
            f'{gc_content*100}%\n'
        )

    return (contig,
            window_start,
            window_stop,
            coverage_short,
            coverage_long,
            gc_content,
            num_frags)





def trim_coverage(window_data:np.ndarray, trim_percentile:int=10):
    """
    function to trim lowest 10% of bins by coverage. If a window is
    below the 10th percentile, coverages and gc are set to NaN and
    num_frags is set to 0
    """
    ten_percentile = np.percentile(window_data['num_frags'], trim_percentile)
    trimmed = window_data.copy()
    in_percentile = window_data['num_frags']<ten_percentile
    trimmed['short'][in_percentile] = np.NaN
    trimmed['long'][in_percentile] = np.NaN
    trimmed['gc'][in_percentile] = np.NaN
    trimmed['num_frags'][in_percentile] = 0
    return trimmed


def delfi(input_file: str,  # TODO: allow AlignmentFile to be used
          autosomes: str,
          reference_file: str,
          blacklist_file: str=None,
          tcmere_file: str=None,
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
    tcmere_file: str
        Path string to a BED4+ file where each interval is a centromere or
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
        stderr.write(f"""
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
        stderr.write(f'Reading genome file...\n')

    contigs = []
    with open(autosomes) as genome:
        for line in genome:
            contents = line.split('\t')
            # account for empty lines
            if len(contents) > 1:
                contigs.append((contents[0],  int(contents[1])))

    tcmeres = []
    if (tcmere_file is not None):
        # TODO: find a standard way to get centromeres, like a track on
        # UCSC
        with open(tcmere_file) as tcmere_bed:
            for line in tcmere_bed:
                region_contig, region_start, region_stop, name, *_ = line.split()
                region_start = int(region_start)
                region_stop = int(region_stop)
                if (contig == region_contig
                    and (name == 'centromere'
                         or name == 'telomere')
                ):
                    tcmeres.append((region_contig, region_start,region_stop, name))

    if verbose:
        stderr.write(f'Generating windows\n')

    # generate DELFI windows
    window_args = []
    for contig, size in contigs:
        for coordinate in range(0, size, window_size):
            # (contig, start, stop)
            window_args.append((
                input_file,
                reference_file,
                contig,
                coordinate,
                coordinate + window_size,
                blacklist_file,
                tcmeres,
                quality_threshold,
                verbose - 1 if verbose > 1 else 0))

    if (verbose):
        stderr.write(f'{len(window_args)} windows created.\n')
        stderr.write('Calculating fragment lengths...\n')

    # pool process to find window frag coverages, gc content
    with Pool(workers) as pool:
        windows = pool.starmap(_delfi_single_window, window_args)

    # move to structured array
    window_array = np.array(
        windows,
        dtype=[('contig', '<U32'),
           ('start', 'u8'),
           ('stop', 'u8'),
           ('short', 'f8'),
           ('long', 'f8'),
           ('gc', 'f8'),
           ('num_frags', 'u8')
        ]
    )

    # remove bottom 10 percentile
    trimmed_windows = trim_coverage(window_array, 10)

    # output
    def _write_out(out: TextIO):
        out.write('#contig\tstart\tstop\tshort\tlong\tgc%\tfrag_count\n')
        for window in trimmed_windows:
            out.write(
                f'{window[0]}\t{window[1]}\t{window[2]}\t{window[3]}\t'
                f'{window[4]}\t{window[5]}\t{window[6]}\n')

    if output_file.endswith('.tsv'):
        with open(output_file, 'w') as out:
            out.write('contig\tstart\tstop\tshort\tlong\tgc%\tfrag_count\n')
            for window in trimmed_windows:
                out.write(
                    f'{window[0]}\t{window[1]}\t{window[2]}\t{window[3]}\t'
                    f'{window[4]}\t{window[5]}\t{window[6]}\n')
    elif output_file.endswith('.bed'):
        with open(output_file, 'w') as out:
            _write_out(out)
    elif output_file.endswith('.bed.gz'):
        with gzip.open(output_file, 'w') as out:
            _write_out(out)
    elif output_file == '-':
        with stdout as out:
            _write_out(out)
    else:
        raise ValueError(
            'Invalid file type! Only .bed, .bed.gz, and .tsv suffixes allowed.'
        )

    num_frags = sum(window[3] for window in windows)

    if (verbose):
        end_time = time.time()
        stderr.write(f'{num_frags} fragments included.\n')
        stderr.write(f'delfi took {end_time - start_time} s to complete\n')
    return None
