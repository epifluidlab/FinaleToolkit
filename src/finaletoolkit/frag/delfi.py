from __future__ import annotations
import time
from multiprocessing.pool import Pool
from typing import Union, TextIO
from sys import stderr, stdout
import gzip

import py2bit
import numpy as np
import pandas
from tqdm import tqdm

from finaletoolkit.frag.delfi_gc_correct import delfi_gc_correct
from finaletoolkit.frag.delfi_merge_bins import delfi_merge_bins
from finaletoolkit.utils.utils import frag_generator, overlaps
from finaletoolkit.genome.gaps import GenomeGaps, ContigGaps


def _delfi_single_window(
        input_file: str,
        reference_file: str,
        contig_gaps: ContigGaps,
        window_start: int,
        window_stop: int,
        blacklist_file: str=None,
        quality_threshold: int=30,
        verbose: Union[int,bool]=False) -> tuple:
    """
    Calculates short and long counts for one window.
    """
    contig = contig_gaps.contig

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

    # check if interval in centromere or telomere
    in_tcmere = contig_gaps.in_tcmere(window_start, window_stop)
    if in_tcmere:
        return (contig,
            window_start,
            window_stop,
            'NOARM',
            np.NaN,
            np.NaN,
            np.NaN,
            0)

    arm = contig_gaps.get_arm(window_start, window_stop)
    # if in short arm
    if arm == 'NOARM':
        return (contig,
            window_start,
            window_stop,
            'NOARM',
            np.NaN,
            np.NaN,
            np.NaN,
            0)
    
    # Iterating on each read in file in specified contig/chromosome
    for _, frag_start, frag_stop, _, _ in frag_generator(
        input_file,
        contig,
        quality_threshold,
        window_start,
        window_stop,
        fraction_low=100,
        fraction_high=220):

        frag_length = frag_stop - frag_start

        assert frag_length > 0, (f"Frag length of {frag_length} found at"
            f"{contig}:{frag_start}-{frag_stop}.")

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
        in_tcmere = contig_gaps.in_tcmere(frag_start, frag_stop)
        if in_tcmere:
            continue

        if (not blacklisted
            and not in_tcmere
        ):
            # append length of fragment to list
            if (frag_length >= 151):
                long_lengths.append(abs(frag_length))
            else:
                short_lengths.append(abs(frag_length))

            frag_pos.append((frag_start, frag_stop))

            num_frags += 1

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
            arm,
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


def delfi(input_file: str,
          autosomes: str,
          bins_file: str,
          reference_file: str,
          blacklist_file: str=None,
          gap_file: Union(str, GenomeGaps)=None,
          output_file: str=None,
          gc_correct:bool=True,
          merge_bins:bool=True,
          window_size: int=5000000,
          subsample_coverage: float=2,
          quality_threshold: int=30,
          workers: int=1,
          preprocessing: bool=True,
          verbose: Union[int, bool]=False) -> pandas.DataFrame:
    """
    A function that replicates the methodology of Christiano et al
    (2019).

    Parameters
    ----------
    input_file: str
        Path string pointing to a bam file containing PE
        fragment reads.
    autosomes: str
        Path string to a .genome file containing only autosomal
        chromosomes
    bins_file: str
        Path string to a BED file containing 100kb bins for reference
        genome of choice. Cristiano et al uses 
    reference_file: str
        Path string to .2bit file.
    blacklist_file: str
        Path string to BED file containing genome blacklist.
    gap_file: str
        Path string to a BED4+ file where each interval is a centromere
        or telomere. A bed file can be used **only if** the fourth field
        for each entry corresponding to a telomere or centromere is
        labled "telomere" or "centromere, respectively.
    output_file: str, optional
        Path to output tsv.
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
        gc_correct: {gc_correct}
        merge_bins: {merge_bins}
        quality_threshold: {quality_threshold}
        workers: {workers}
        preprocessing: {preprocessing}
        verbose: {verbose}
        \n""")

    if verbose:
        stderr.write(f'Reading genome file...\n')

    # Read chromosome names and lengths from .genome file
    contigs = []
    with open(autosomes) as genome:
        for line in genome:
            contents = line.split('\t')
            # account for empty lines
            if len(contents) > 1:
                contigs.append((contents[0],  int(contents[1])))

    # Prepare genome gaps using GenomeGaps class
    gaps = None
    if (gap_file is not None):
        if type(gap_file) == str:
            gaps = GenomeGaps(gap_file)
        elif type(gap_file) == GenomeGaps:
            gaps = gap_file
        else:
            raise TypeError(
                f'{type(gap_file)} is not accepted type for gap_file'
            )

    # Read 100kb bins and filter out bins that overlap gaps, darkregions
    
    # opening 100kb bins BED file into a dataframe
    if verbose:
        stderr.write(f'Opening bins file...\n')

    bins = pandas.read_csv(
        bins_file,
        names=["contig", "start", "stop"],
        usecols=[0, 1, 2],
        dtype={"contig":str, "start":np.int32, "stop":np.int32},
        delimiter='\t',
        comment='#',
    )

    if verbose:
        stderr.write(f'{bins.shape[0]} bins read from file.\n')
        stderr.write(f'Filtering gaps...\n')
        
    # filtering for gaps
    if gaps is not None:
        # finding overlap
        overlaps_gap = overlaps(
            bins['contig'].to_numpy(),
            bins['start'].to_numpy(),
            bins['stop'].to_numpy(),
            gaps.gaps['contig'],
            gaps.gaps['start'],
            gaps.gaps['stop'],
        )
        # masking by overlap
        gapless_bins = bins.loc[~overlaps_gap]
        if verbose:
            stderr.write(f'{bins.shape[0]-gapless_bins.shape[0]} bins removed'
                         '\n')
    else:
        if verbose:
            stderr.write(f'No gaps specified, skipping.\n')
        gapless_bins = bins

    # filtering for darkregions
    """
    if verbose:
        stderr.write(f'Filtering darkregions...\n')
    if blacklist_file is not None:
        # opening blacklist file into a dataframe
        darkregions = pandas.read_csv(
            blacklist_file,
            names=["contig", "start", "stop"],
            usecols=[0, 1, 2],
            dtype={"contig":str, "start":np.int32, "stop":np.int32},
            delimiter='\t'
        )
        # finding overlap
        overlaps_darkregion = overlaps(
            gapless_bins['contig'].to_numpy(),
            gapless_bins['start'].to_numpy(),
            gapless_bins['stop'].to_numpy(),
            darkregions['contig'].to_numpy(),
            darkregions['start'].to_numpy(),
            darkregions['stop'].to_numpy(),
        )
        # masking by overlap
        darkless_bins = gapless_bins.loc[~overlaps_darkregion]
        if verbose:
            stderr.write(f'{gapless_bins.shape[0]-darkless_bins.shape[0]} bins'
                         ' removed\n')
    else:
        if verbose:
            stderr.write(f'No darkregions given, skipping.\n')
        darkless_bins = gapless_bins
    """
    darkless_bins = gapless_bins

    # generating args for pooled processes
    if verbose:
        stderr.write(f'Preparing to generate short and long coverages.\n')

    window_args = []
    contig_gaps = None

    for contig, size in contigs:
        if gaps is not None:
            # print(contig)
            contig_gaps = gaps.get_contig_gaps(contig)
        else:
            contig_gaps = None
        for _, start, stop, *_ in (
            darkless_bins.loc[darkless_bins.loc[:,'contig']==contig]
            .itertuples(index=False, name=None)
        ):
            window_args.append((
                input_file,
                reference_file,
                contig_gaps,
                start,
                stop,
                blacklist_file,
                quality_threshold,
                verbose - 1 if verbose > 1 else 0))

    if (verbose):
        stderr.write(f'{len(window_args)} windows created.\n')
        stderr.write('Calculating fragment lengths...\n')

    # pool process to find window frag coverages, gc content
    with Pool(workers) as pool:
        windows = pool.starmap(
            _delfi_single_window, tqdm(window_args), 50)


    # move to dataframe
    if (verbose):
        stderr.write('Done.\n')
        stderr.write('Removing remaining accrocentric bins...\n')
   
    window_df = pandas.DataFrame(
        windows,
        columns=[
            'contig', 'start', 'stop', 'arm', 'short', 'long', 'gc',
            'num_frags']
    )
    # remove remaining NOARM bins
    trimmed_windows = window_df.loc[window_df['arm']!='NOARM', :].copy()

    if (verbose):
        stderr.write(f'{trimmed_windows.shape[0]} bins remaining...\n')

    # calculating ratio
    if (verbose):
        stderr.write('Calculating ratio...\n')

    trimmed_windows['ratio'] = trimmed_windows['short']/trimmed_windows['long']

    # gc correct
    if gc_correct:
        if (verbose):
            stderr.write('GC bias correction...\n')
        gc_corrected = delfi_gc_correct(trimmed_windows, 0.75, 8, verbose)
    else:
        gc_corrected = trimmed_windows

    if merge_bins:
        if (verbose):
            stderr.write('Merging bins...\n')
        final_bins = delfi_merge_bins(
            gc_corrected, gc_correct, True, verbose=verbose)
    else:
        final_bins = gc_corrected
    
    # output
    if (verbose):
        stderr.write(f'{final_bins.shape[0]} bins remaining.\n')
        
    output_delfi = final_bins.rename(columns={'contig':'#contig'})

    if output_file.endswith('.bed') or output_file.endswith('.tsv'):
        output_delfi.to_csv(output_file, sep='\t', index=False)
    elif output_file.endswith('.csv'):
        output_delfi.to_csv(output_file, sep=',', index=False)
    elif output_file.endswith('.bed.gz'):
        output_delfi.to_csv(
            output_file,
            sep='\t',
            index=False,
            encoding='gzip')
    elif output_file == '-':
        with stdout as out:
            for window in gc_corrected.itertuples():
                tab_separated = "\t".join(window)
                out.write(
                    f'{tab_separated}\n')
    else:
        raise ValueError(
            'Invalid file type! Only .bed, .bed.gz, and .tsv suffixes allowed.'
        )

    num_frags = sum(window[7] for window in windows)

    if (verbose):
        end_time = time.time()
        stderr.write(f'{num_frags} fragments included.\n')
        stderr.write(f'delfi took {end_time - start_time} s to complete\n')
    return gc_corrected
