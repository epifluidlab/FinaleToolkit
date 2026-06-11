from __future__ import annotations

import time
import warnings
from collections import defaultdict
from multiprocessing.pool import Pool
from sys import stderr, stdout
from typing import Union

import numpy as np
import pandas
from tqdm import tqdm

from finaletoolkit.frag._delfi_gc_correct import delfi_gc_correct
from finaletoolkit.frag._delfi_merge_bins import delfi_merge_bins
from finaletoolkit.genome.gaps import GenomeGaps
from finaletoolkit.io.alignment import AlignmentWrapper
from finaletoolkit.io.reference import ReferenceWrapper
from finaletoolkit.utils.logging import get_logger
from finaletoolkit.utils import chrom_sizes_to_list, overlaps
from finaletoolkit.utils.validation import valid_interval

logger = get_logger(__name__)


def trim_coverage(window_data:np.ndarray, trim_percentile:int=10):
    """
    function to trim lowest 10% of bins by coverage. If a window is
    below the 10th percentile, coverages and gc are set to NaN and
    num_frags is set to 0
    """
    ten_percentile = np.percentile(window_data['num_frags'], trim_percentile)
    trimmed = window_data.copy()
    in_percentile = window_data['num_frags']<ten_percentile
    trimmed['short'][in_percentile] = np.nan
    trimmed['long'][in_percentile] = np.nan
    trimmed['gc'][in_percentile] = np.nan
    trimmed['num_frags'][in_percentile] = 0
    return trimmed


# ---------------------------------------------------------------------------
# Worker-process state. Set once per worker via the Pool initializer so that
# files are opened and inputs parsed once per worker instead of once per
# window. The previous implementation reopened the alignment and reference
# files and re-parsed the blacklist BED for every 100kb window (~26K windows
# per whole genome), which dominated runtime.
# ---------------------------------------------------------------------------
_WORKER_ALIGNMENT = None  # AlignmentWrapper (BAM/CRAM/frag.gz)
_WORKER_REF = None        # ReferenceWrapper (.2bit/FASTA)
_WORKER_BLACKLIST = None  # dict[contig] -> (sorted_starts, stops_aligned_to_starts)
_WORKER_CONTIG_GAPS = None  # dict[contig] -> ContigGaps


def _delfi_pool_initializer(input_file, reference_file, quality_threshold,
                            blacklist_by_contig, contig_gaps_by_contig):
    """Open shared handles and stash preparsed inputs in worker globals."""
    global _WORKER_ALIGNMENT, _WORKER_REF
    global _WORKER_BLACKLIST, _WORKER_CONTIG_GAPS
    # AlignmentWrapper transparently handles BAM, CRAM, and tabix-indexed
    # frag.gz/bed.gz input, mirroring frag_generator's input handling.
    _WORKER_ALIGNMENT = AlignmentWrapper(
        input_file,
        reference_file=reference_file,
        quality_threshold=quality_threshold,
    )
    # Each worker is single-threaded, so the reference does not need a lock.
    _WORKER_REF = ReferenceWrapper(reference_file, use_lock=False)
    _WORKER_BLACKLIST = blacklist_by_contig
    _WORKER_CONTIG_GAPS = contig_gaps_by_contig


def _load_blacklist_indexed(blacklist_file):
    """Parse the blacklist BED once and index it by contig.

    Returns a dict mapping each contig to a pair of sorted, position-aligned
    numpy arrays ``(starts, stops)``. This is built once in the parent process
    and shared with every worker, replacing the per-window re-read and linear
    scan of the original implementation.
    """
    if blacklist_file is None:
        return {}
    by_contig = defaultdict(list)
    with open(blacklist_file) as blacklist:
        for line in blacklist:
            parts = line.split()
            if len(parts) < 3:
                continue
            by_contig[parts[0]].append((int(parts[1]), int(parts[2])))
    out = {}
    for contig, regions in by_contig.items():
        regions.sort()
        starts = np.array([r[0] for r in regions], dtype=np.int64)
        stops = np.array([r[1] for r in regions], dtype=np.int64)
        out[contig] = (starts, stops)
    return out


def _blacklist_in_window(contig, window_start, window_stop):
    """Return blacklist regions fully contained in ``[window_start, window_stop]``.

    Equivalent to the original ``window_start <= region_start and
    window_stop >= region_stop`` filter, but uses binary search over the
    contig's sorted start positions instead of scanning the whole file.
    """
    if not _WORKER_BLACKLIST or contig not in _WORKER_BLACKLIST:
        return ()
    starts, stops = _WORKER_BLACKLIST[contig]
    lo = np.searchsorted(starts, window_start, side='left')
    if lo >= len(starts):
        return ()
    sub_starts = starts[lo:]
    sub_stops = stops[lo:]
    keep = sub_stops <= window_stop
    return tuple(zip(sub_starts[keep].tolist(), sub_stops[keep].tolist()))


def delfi(input_file: str,
          chrom_sizes: str,
          bins_file: str,
          reference_file: str,
          blacklist_file: str=None,
          gap_file: Union[str, GenomeGaps]=None,
          output_file: str=None,
          no_gc_correct: bool=False,
          gc_correct: bool | None = None,
          remove_nocov:bool=True,
          merge_bins:bool=True,
          window_size: int=5000000,
          quality_threshold: int=30,
          workers: int=1,
          verbose: Union[int, bool]=False) -> pandas.DataFrame:
    """
    A function that replicates the methodology of Christiano et al
    (2019).

    Parameters
    ----------
    input_file: str
        Path string pointing to a BAM, CRAM, or tabix-indexed fragment
        (.frag.gz/.bed.gz) file containing PE fragment reads.
    chrom_sizes: str
        Path string to a chrom.sizes file containing only autosomal
        chromosomes
    bins_file: str
        Path string to a BED file containing 100kb bins for reference
        genome of choice.
    reference_file: str
        Path string to reference genome file (.2bit or FASTA). Used for
        GC-content correction. When `input_file` is a CRAM file, this
        must be a FASTA file (not .2bit) because htslib requires FASTA
        for CRAM decoding.
    gap_file: str or GenomeGaps
        Specifies locations of telomeres and centromeres for reference 
        genome. There are three options:
        - Path string to a BED4+ file where each interval is a
        centromere or telomere. A bed file can be used **only if** the 
        fourth field for each entry corresponding to a telomere or
        centromere is labled "telomere" or "centromere, respectively. 
        - String naming reference genome used. Options are "b37",
        "hg19", "hg38", and "GRCh38".
        - Alternatively, a finaletoolkit.genome.GenomeGaps with gap 
        info associated with the reference genome of choice may be
        used.
    blacklist_file: str
        Path string to BED file containing genome blacklist regions.
    output_file: str, optional
        Path to output tsv.
    no_gc_correct: bool
        Skip gc-correction. Default is False.
    gc_correct: bool, optional
        Deprecated command to perform gc-correction. Use no_gc_correct instead. 
    remove_nocov: bool
        Remove two windows described by Cristiano et al (2019) as low
        coverage. These windows might not apply to reference genomes
        other than hg19. Default is True.
    merge_bins: bool
        Perform merging from 100kb bins to 5Mb bins. Default is True.
    window_size: int
        Size (in bases) of non-overlapping windows to cover genome. Default is
        5000000.
    workers: int, optional
        Number of worker processes to use. Default is 1.
    verbose: int or bool, optional
        Determines how many print statements and loading bars appear in
        stdout. Default is False.

    Returns
    -------
    pandas DataFrame
        Results of delfi analysis, with column names corresponding to
        those generated by the original author's scripts.
    """

    if (verbose):
        start_time = time.time()
        stderr.write(f"""
        input_file: {input_file}
        chrom_sizes: {chrom_sizes}
        reference_file: {reference_file}
        blacklist_file: {blacklist_file}
        output_file: {output_file}
        window_size: {window_size}
        no_gc_correct: {no_gc_correct}
        gc_correct: {gc_correct} (overrides no_gc_correct if set)
        remove_nocov: {remove_nocov}
        merge_bins: {merge_bins}
        quality_threshold: {quality_threshold}
        workers: {workers}
        verbose: {verbose}
        \n""")

    if verbose:
        stderr.write('Reading genome file...\n')

    # Read chromosome names and lengths from .genome file
    contigs = chrom_sizes_to_list(chrom_sizes)

    # Flip boolean value of no_gc_correct to use gc_correct in following logic if not set
    # If set, use the value but give a deprecation warning
    if gc_correct is None:
        gc_correct = not no_gc_correct
    else:
        warnings.warn("Warning: gc_correct is deprecated and may be removed in future releases. Use no_gc_correct instead")

    # Prepare genome gaps using GenomeGaps class
    gaps = None
    if isinstance(gap_file, str):
        gaps = GenomeGaps(gap_file)
    elif isinstance(gap_file, GenomeGaps):
        gaps = gap_file
    elif gaps is None:
        pass
    else:
        raise TypeError(
            f'{type(gap_file)} is not accepted type for gap_file'
        )

    # Read 100kb bins and filter out bins that overlap gaps, darkregions

    # opening 100kb bins BED file into a dataframe
    if verbose:
        stderr.write('Opening bins file...\n')

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
        stderr.write('Filtering gaps...\n')

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
            stderr.write(f'{bins.shape[0]-gapless_bins.shape[0]} bins '
            'removed\n')
    else:
        if verbose:
            stderr.write('No gaps specified, skipping.\n')
        gapless_bins = bins

    # generating args for pooled processes
    if verbose:
        stderr.write('Preparing to generate short and long coverages.\n')

    # Pre-parse shared, read-only inputs once. These are handed to each
    # worker via the Pool initializer rather than pickled into every task.
    blacklist_by_contig = _load_blacklist_indexed(blacklist_file)
    contig_gaps_by_contig = {}
    if gaps is not None:
        for contig, _size in contigs:
            contig_gaps_by_contig[contig] = gaps.get_contig_gaps(contig)

    window_args = []
    for contig, _size in contigs:
        for _, start, stop, *_ in (
            gapless_bins.loc[gapless_bins.loc[:,'contig']==contig]
            .itertuples(index=False, name=None)
        ):
            window_args.append((
                contig,
                start,
                stop,
                verbose - 1 if verbose > 1 else 0))

    if (verbose):
        stderr.write(f'{len(window_args)} windows created.\n')
        stderr.write('Calculating fragment lengths...\n')

    # pool process to find window frag coverages, gc content
    with Pool(
        workers,
        initializer=_delfi_pool_initializer,
        initargs=(input_file, reference_file, quality_threshold,
                  blacklist_by_contig, contig_gaps_by_contig),
    ) as pool:
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

    trimmed_windows['ratio'] = np.where(trimmed_windows['long'] == 0, np.nan, trimmed_windows['short'] / trimmed_windows['long']) #handle the case where the 'long' is 0.

    # remove nocov windows
    if remove_nocov:
        no_nocov_slice = np.logical_and(
            np.logical_not(np.equal(np.arange(trimmed_windows.shape[0]),8779)),
            np.logical_not(np.equal(np.arange(trimmed_windows.shape[0]), 13664)))
        corrected_delfi_drop_nocov = trimmed_windows.loc[no_nocov_slice].reset_index()
    else:
        corrected_delfi_drop_nocov = trimmed_windows

    # gc correct
    if gc_correct:
        if (verbose):
            stderr.write('GC bias correction...\n')
        gc_corrected = delfi_gc_correct(corrected_delfi_drop_nocov, 0.75, 8, verbose)
    else:
        gc_corrected = corrected_delfi_drop_nocov

    # merge bins
    if merge_bins:
        if (verbose):
            stderr.write('Merging bins...\n')
        final_bins = delfi_merge_bins(
            gc_corrected, gc_correct, verbose=verbose)
    else:
        final_bins = gc_corrected
    
    # output
    if (verbose):
        stderr.write(f'{final_bins.shape[0]} bins remaining.\n')

    if output_file is not None:
        output_delfi = final_bins.rename(columns={'contig':'#contig'})
        if output_file.endswith('.bed') or output_file.endswith('.tsv'):
            output_delfi.to_csv(output_file, sep='\t', index=False)
        elif output_file.endswith('.csv'):
            final_bins.to_csv(output_file, sep=',', index=False)
        elif output_file.endswith('.bed.gz'):
            output_delfi.to_csv(
                output_file,
                sep='\t',
                index=False,
                encoding='gzip')
        elif output_file == '-':
            with stdout as out:
                for window in final_bins.itertuples():
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
    return final_bins


def _delfi_single_window(
        contig: str,
        window_start: int,
        window_stop: int,
        verbose: Union[int,bool]=False) -> tuple:
    """
    Calculates short and long counts for one window.

    Reads the alignment, reference, blacklist, and gap inputs from the
    worker-process globals set by ``_delfi_pool_initializer`` so that no
    file is reopened and no input re-parsed per window.
    """

    contig_gaps = (
        _WORKER_CONTIG_GAPS.get(contig)
        if _WORKER_CONTIG_GAPS is not None else None
    )

    if contig_gaps is not None:
        # check if interval in centromere or telomere
        if contig_gaps.in_tcmere(window_start, window_stop):
            return (contig,
                window_start,
                window_stop,
                'NOARM',
                np.nan,
                np.nan,
                np.nan,
                0)

        arm = contig_gaps.get_arm(window_start, window_stop)
        # if in short arm
        if arm == 'NOARM':
            return (contig,
                window_start,
                window_stop,
                'NOARM',
                np.nan,
                np.nan,
                np.nan,
                0)
    else:
        arm = contig

    blacklist_regions = _blacklist_in_window(contig, window_start, window_stop)

    short_lengths = []
    long_lengths = []
    num_frags = 0

    # Iterating on each fragment in the window. The shared AlignmentWrapper
    # applies the same quality / read-pair filtering as frag_generator; we
    # replicate frag_generator's length window (100-220) and its default
    # "midpoint" intersect policy here so output is unchanged.
    for frag in _WORKER_ALIGNMENT.fetch(contig, window_start, window_stop):
        frag_start = frag.start
        frag_stop = frag.stop
        frag_length = frag_stop - frag_start

        if frag_length < 100 or frag_length > 220:
            continue

        midpoint = (frag_start + frag_stop) // 2
        if midpoint < window_start or midpoint >= window_stop:
            continue

        # check if in blacklist
        blacklisted = False
        for region in blacklist_regions:
            if (
                (frag_start >= region[0] and frag_start < region[1])
                and (frag_stop >= region[0] and frag_stop < region[1])
            ):
                blacklisted = True
                break
        if blacklisted:
            continue

        # check if in centromere or telomere
        if contig_gaps is not None and contig_gaps.in_tcmere(frag_start, frag_stop):
            continue

        # append length of fragment to list
        if (frag_length >= 151):
            long_lengths.append(abs(frag_length))
        else:
            short_lengths.append(abs(frag_length))

        num_frags += 1

    # finding gc amount
    if not valid_interval(
        _WORKER_REF.chroms, contig, window_start, window_stop,
        throw_on_error=False,
    ):
        logger.warning(
            f"Invalid interval {contig}:{window_start}-{window_stop} for "
            "reference. Skipping GC calculation."
        )
        ref_bases = ""
    else:
        ref_bases = _WORKER_REF.sequence(contig, window_start, window_stop)

    # ref_bases is upper-cased by ReferenceWrapper, so counting 'G'/'C'
    # directly matches the previous per-base comparison.
    num_gc = ref_bases.count('G') + ref_bases.count('C')

    # window_length
    window_coverage = window_stop - window_start

    # NaN if no fragments in window.
    gc_content = num_gc / window_coverage if num_frags > 0 else np.nan

    coverage_short = len(short_lengths)
    coverage_long = len(long_lengths)

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
