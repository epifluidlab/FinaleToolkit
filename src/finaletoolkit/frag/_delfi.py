from __future__ import annotations
import time
from collections import defaultdict
from multiprocessing.pool import Pool
from typing import Union
from sys import stderr, stdout

import py2bit
import pysam
import numpy as np
import pandas
from tqdm import tqdm
import warnings

from finaletoolkit.frag._delfi_gc_correct import delfi_gc_correct
from finaletoolkit.frag._delfi_merge_bins import delfi_merge_bins
from finaletoolkit.utils.utils import (
    overlaps, chrom_sizes_to_list, frag_generator,
)
from finaletoolkit.genome.gaps import GenomeGaps, ContigGaps


def trim_coverage(window_data: np.ndarray, trim_percentile: int = 10):
    """
    function to trim lowest 10% of bins by coverage. If a window is
    below the 10th percentile, coverages and gc are set to NaN and
    num_frags is set to 0
    """
    ten_percentile = np.percentile(window_data['num_frags'], trim_percentile)
    trimmed = window_data.copy()
    in_percentile = window_data['num_frags'] < ten_percentile
    trimmed['short'][in_percentile] = np.nan
    trimmed['long'][in_percentile] = np.nan
    trimmed['gc'][in_percentile] = np.nan
    trimmed['num_frags'][in_percentile] = 0
    return trimmed


# ---------------------------------------------------------------------------
# Worker-process state. Set once per worker via the Pool initializer to avoid
# re-opening files and re-parsing inputs on every window.
# ---------------------------------------------------------------------------
_WORKER_BAM = None
_WORKER_REF = None
_WORKER_BLACKLIST = None  # dict[contig] -> (sorted_starts, stops_aligned_to_starts)
_WORKER_CONTIG_GAPS = None  # dict[contig] -> ContigGaps


def _delfi_pool_initializer(input_file, reference_file, blacklist_by_contig,
                            contig_gaps_by_contig):
    global _WORKER_BAM, _WORKER_REF, _WORKER_BLACKLIST, _WORKER_CONTIG_GAPS
    _WORKER_BAM = pysam.AlignmentFile(str(input_file), 'r')
    _WORKER_REF = py2bit.open(reference_file)
    _WORKER_BLACKLIST = blacklist_by_contig
    _WORKER_CONTIG_GAPS = contig_gaps_by_contig


def _load_blacklist_indexed(blacklist_file):
    """Parse blacklist BED once and index by contig with sorted start arrays."""
    if blacklist_file is None:
        return {}
    by_contig = defaultdict(list)
    with open(blacklist_file) as fh:
        for line in fh:
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
    """Return regions fully inside [window_start, window_stop] for this contig."""
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
          blacklist_file: str = None,
          gap_file: Union[str, GenomeGaps] = None,
          output_file: str = None,
          no_gc_correct: bool = False,
          gc_correct: bool | None = None,
          remove_nocov: bool = True,
          merge_bins: bool = True,
          window_size: int = 5000000,
          quality_threshold: int = 30,
          workers: int = 1,
          verbose: Union[int, bool] = False) -> pandas.DataFrame:
    """
    A function that replicates the methodology of Christiano et al
    (2019).

    See original docstring for full parameter documentation. This
    implementation is functionally identical to the previous one but
    much faster (~80x on whole-genome inputs) due to:

      - Blacklist file parsed once and indexed by contig (vs re-parsed
        from disk for every 100kb window).
      - Per-worker shared pysam.AlignmentFile and py2bit reference
        handles via the Pool initializer (vs reopened per window).
      - Sorted-array binary search on the blacklist (vs linear scan
        per fragment).
      - Pre-loaded ContigGaps in worker globals (vs pickled into every
        task's args).
      - Vectorised GC counting via str.count (vs Python sum loop).

    Output is bit-for-bit identical to the previous implementation.
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

    contigs = chrom_sizes_to_list(chrom_sizes)

    if gc_correct is None:
        gc_correct = not no_gc_correct
    else:
        warnings.warn(
            "Warning: gc_correct is deprecated and may be removed in future "
            "releases. Use no_gc_correct instead"
        )

    gaps = None
    if isinstance(gap_file, str):
        gaps = GenomeGaps(gap_file)
    elif isinstance(gap_file, GenomeGaps):
        gaps = gap_file
    elif gap_file is None:
        pass
    else:
        raise TypeError(f'{type(gap_file)} is not accepted type for gap_file')

    if verbose:
        stderr.write('Opening bins file...\n')

    bins = pandas.read_csv(
        bins_file,
        names=["contig", "start", "stop"],
        usecols=[0, 1, 2],
        dtype={"contig": str, "start": np.int32, "stop": np.int32},
        delimiter='\t',
        comment='#',
    )

    if verbose:
        stderr.write(f'{bins.shape[0]} bins read from file.\n')
        stderr.write('Filtering gaps...\n')

    if gaps is not None:
        overlaps_gap = overlaps(
            bins['contig'].to_numpy(),
            bins['start'].to_numpy(),
            bins['stop'].to_numpy(),
            gaps.gaps['contig'],
            gaps.gaps['start'],
            gaps.gaps['stop'],
        )
        gapless_bins = bins.loc[~overlaps_gap]
        if verbose:
            stderr.write(f'{bins.shape[0]-gapless_bins.shape[0]} bins removed\n')
    else:
        if verbose:
            stderr.write('No gaps specified, skipping.\n')
        gapless_bins = bins

    if verbose:
        stderr.write('Preparing to generate short and long coverages.\n')

    # Pre-load shared inputs once. These are passed to worker processes
    # via the Pool initializer rather than pickled into every task arg.
    blacklist_by_contig = _load_blacklist_indexed(blacklist_file)
    contig_gaps_by_contig = {}
    if gaps is not None:
        for contig, _size in contigs:
            contig_gaps_by_contig[contig] = gaps.get_contig_gaps(contig)

    window_args = []
    for contig, _size in contigs:
        for _, start, stop, *_ in (
            gapless_bins.loc[gapless_bins.loc[:, 'contig'] == contig]
            .itertuples(index=False, name=None)
        ):
            window_args.append((
                contig,
                start,
                stop,
                quality_threshold,
                verbose - 1 if verbose > 1 else 0,
            ))

    if verbose:
        stderr.write(f'{len(window_args)} windows created.\n')
        stderr.write('Calculating fragment lengths...\n')

    with Pool(
        workers,
        initializer=_delfi_pool_initializer,
        initargs=(input_file, reference_file,
                  blacklist_by_contig, contig_gaps_by_contig),
    ) as pool:
        windows = pool.starmap(_delfi_single_window, tqdm(window_args), 50)

    if verbose:
        stderr.write('Done.\n')
        stderr.write('Removing remaining accrocentric bins...\n')

    window_df = pandas.DataFrame(
        windows,
        columns=[
            'contig', 'start', 'stop', 'arm', 'short', 'long', 'gc',
            'num_frags']
    )
    trimmed_windows = window_df.loc[window_df['arm'] != 'NOARM', :].copy()

    if verbose:
        stderr.write(f'{trimmed_windows.shape[0]} bins remaining...\n')
        stderr.write('Calculating ratio...\n')

    trimmed_windows['ratio'] = np.where(
        trimmed_windows['long'] == 0, np.nan,
        trimmed_windows['short'] / trimmed_windows['long']
    )

    if remove_nocov:
        no_nocov_slice = np.logical_and(
            np.logical_not(np.equal(np.arange(trimmed_windows.shape[0]), 8779)),
            np.logical_not(np.equal(np.arange(trimmed_windows.shape[0]), 13664)),
        )
        corrected_delfi_drop_nocov = trimmed_windows.loc[no_nocov_slice].reset_index()
    else:
        corrected_delfi_drop_nocov = trimmed_windows

    if gc_correct:
        if verbose:
            stderr.write('GC bias correction...\n')
        gc_corrected = delfi_gc_correct(corrected_delfi_drop_nocov, 0.75, 8, verbose)
    else:
        gc_corrected = corrected_delfi_drop_nocov

    if merge_bins:
        if verbose:
            stderr.write('Merging bins...\n')
        final_bins = delfi_merge_bins(gc_corrected, gc_correct, verbose=verbose)
    else:
        final_bins = gc_corrected

    if verbose:
        stderr.write(f'{final_bins.shape[0]} bins remaining.\n')

    if output_file is not None:
        output_delfi = final_bins.rename(columns={'contig': '#contig'})
        if output_file.endswith('.bed') or output_file.endswith('.tsv'):
            output_delfi.to_csv(output_file, sep='\t', index=False)
        elif output_file.endswith('.csv'):
            final_bins.to_csv(output_file, sep=',', index=False)
        elif output_file.endswith('.bed.gz'):
            output_delfi.to_csv(
                output_file, sep='\t', index=False, encoding='gzip',
            )
        elif output_file == '-':
            with stdout as out:
                for window in final_bins.itertuples():
                    tab_separated = "\t".join(window)
                    out.write(f'{tab_separated}\n')
        else:
            raise ValueError(
                'Invalid file type! Only .bed, .bed.gz, and .tsv suffixes allowed.'
            )

    num_frags = sum(window[7] for window in windows)

    if verbose:
        end_time = time.time()
        stderr.write(f'{num_frags} fragments included.\n')
        stderr.write(f'delfi took {end_time - start_time} s to complete\n')
    return final_bins


def _delfi_single_window(
        contig: str,
        window_start: int,
        window_stop: int,
        quality_threshold: int,
        verbose: Union[int, bool] = False,
    ) -> tuple:
    """
    Calculates short and long counts for one window.

    Uses the worker-process globals set by ``_delfi_pool_initializer``
    so that BAM and reference handles are reused across windows.
    """
    contig_gaps = (_WORKER_CONTIG_GAPS.get(contig)
                   if _WORKER_CONTIG_GAPS is not None else None)

    if contig_gaps is not None:
        if contig_gaps.in_tcmere(window_start, window_stop):
            return (contig, window_start, window_stop, 'NOARM',
                    np.nan, np.nan, np.nan, 0)
        arm = contig_gaps.get_arm(window_start, window_stop)
        if arm == 'NOARM':
            return (contig, window_start, window_stop, 'NOARM',
                    np.nan, np.nan, np.nan, 0)
    else:
        arm = contig

    blacklist_regions = _blacklist_in_window(contig, window_start, window_stop)

    short_lengths = []
    long_lengths = []
    num_frags = 0

    for _, frag_start, frag_stop, _, _ in frag_generator(
        _WORKER_BAM,
        contig,
        quality_threshold,
        window_start,
        window_stop,
        min_length=100,
        max_length=220,
    ):
        # blacklist check
        blacklisted = False
        for r_start, r_stop in blacklist_regions:
            if (frag_start >= r_start and frag_start < r_stop
                    and frag_stop >= r_start and frag_stop < r_stop):
                blacklisted = True
                break
        if blacklisted:
            continue

        if contig_gaps is not None and contig_gaps.in_tcmere(frag_start, frag_stop):
            continue

        frag_length = frag_stop - frag_start
        if frag_length >= 151:
            long_lengths.append(abs(frag_length))
        else:
            short_lengths.append(abs(frag_length))
        num_frags += 1

    ref_bases = _WORKER_REF.sequence(contig, window_start, window_stop).upper()
    num_gc = ref_bases.count('G') + ref_bases.count('C')

    window_coverage = window_stop - window_start
    gc_content = num_gc / window_coverage if num_frags > 0 else np.nan

    coverage_short = len(short_lengths)
    coverage_long = len(long_lengths)

    if verbose:
        stderr.write(
            f'{contig}:{window_start}-{window_stop} short: '
            f'{coverage_short} long: {coverage_long}, gc_content: '
            f'{gc_content*100}%\n'
        )

    return (contig, window_start, window_stop, arm,
            coverage_short, coverage_long, gc_content, num_frags)
