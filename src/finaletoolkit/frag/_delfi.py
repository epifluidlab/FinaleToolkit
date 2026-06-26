"""
DELFI fragmentation feature (Cristiano et al., 2019): short/long fragment
ratios per genomic bin, with optional GC correction and 100kb->5Mb merging.
"""
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
from finaletoolkit.utils import chrom_sizes_to_list, overlaps
from finaletoolkit.utils.logging import get_logger
from finaletoolkit.utils.validation import valid_interval

logger = get_logger(__name__)

__all__ = ["delfi", "trim_coverage"]


def trim_coverage(window_data: np.ndarray, trim_percentile: int = 10):
    """Trim the lowest ``trim_percentile``% of bins by fragment count.

    Bins below the percentile have ``short``/``long``/``gc`` set to NaN and
    ``num_frags`` set to 0.
    """
    threshold = np.percentile(window_data["num_frags"], trim_percentile)
    trimmed = window_data.copy()
    in_percentile = window_data["num_frags"] < threshold
    trimmed["short"][in_percentile] = np.nan
    trimmed["long"][in_percentile] = np.nan
    trimmed["gc"][in_percentile] = np.nan
    trimmed["num_frags"][in_percentile] = 0
    return trimmed


# ---------------------------------------------------------------------------
# Worker-process state, set once per worker via the Pool initializer so files
# are opened and inputs parsed once per worker instead of once per window. The
# previous implementation reopened the alignment and reference files and
# re-parsed the blacklist BED for every 100kb window (~26K windows per whole
# genome), which dominated runtime.
#
# This worker-pool optimization (~80x faster, identical output) was
# contributed by D.H.K. (Duco) Gaillard (@DucoG) in upstream PR #172
# (epifluidlab/FinaleToolkit#172) and ported onto this implementation.
# ---------------------------------------------------------------------------
_WORKER_ALIGNMENT = None  # AlignmentWrapper (BAM/CRAM/frag.gz)
_WORKER_REF = None  # ReferenceWrapper (.2bit/FASTA)
_WORKER_BLACKLIST = None  # dict[contig] -> (sorted_starts, stops)
_WORKER_CONTIG_GAPS = None  # dict[contig] -> ContigGaps


def _delfi_pool_initializer(
    input_file, reference_file, quality_threshold,
    blacklist_by_contig, contig_gaps_by_contig,
):
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
    # Each worker is single-threaded, so the reference needs no lock.
    _WORKER_REF = ReferenceWrapper(reference_file, use_lock=False)
    _WORKER_BLACKLIST = blacklist_by_contig
    _WORKER_CONTIG_GAPS = contig_gaps_by_contig


def _load_blacklist_indexed(blacklist_file):
    """Parse the blacklist BED once and index it by contig.

    Returns a dict mapping each contig to a pair of sorted, position-aligned
    numpy arrays ``(starts, stops)``. Built once in the parent process and
    shared with every worker, replacing the per-window re-read and linear scan.
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
    """Blacklist regions fully contained in ``[window_start, window_stop]``.

    Equivalent to the original ``window_start <= region_start and
    window_stop >= region_stop`` filter, but uses binary search over the
    contig's sorted start positions instead of scanning the whole file.
    """
    if not _WORKER_BLACKLIST or contig not in _WORKER_BLACKLIST:
        return ()
    starts, stops = _WORKER_BLACKLIST[contig]
    lo = np.searchsorted(starts, window_start, side="left")
    if lo >= len(starts):
        return ()
    sub_starts = starts[lo:]
    sub_stops = stops[lo:]
    keep = sub_stops <= window_stop
    return tuple(zip(sub_starts[keep].tolist(), sub_stops[keep].tolist()))


def delfi(
    input_file: str,
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
    verbose: Union[int, bool] = False,
) -> pandas.DataFrame:
    """Compute DELFI features, replicating Cristiano et al. (2019).

    Parameters
    ----------
    input_file : str
        BAM/CRAM/fragment file of paired-end reads.
    chrom_sizes : str
        ``.chrom.sizes`` file (autosomes only to match the original scripts).
    bins_file : str
        BED of bins (100kb to match the original methodology).
    reference_file : str
        Reference genome (``.2bit`` or FASTA; FASTA required for CRAM input).
    blacklist_file : str, optional
        BED of regions to ignore.
    gap_file : str or GenomeGaps, optional
        Telomere/centromere annotations: a BED4 file, a genome name
        (``"b37"``/``"hg19"``/``"hg38"``/``"GRCh38"``), or a
        :class:`~finaletoolkit.genome.GenomeGaps`.
    output_file : str, optional
        Output BED/`.bed.gz`/`.tsv`/`.csv` path, or ``"-"`` for stdout.
    no_gc_correct : bool, optional
        Skip GC correction (default ``False``).
    gc_correct : bool, optional
        Deprecated; overrides ``no_gc_correct`` when set.
    remove_nocov : bool, optional
        Remove the two hg19 no-coverage windows of Cristiano et al. (default
        ``True``).  Disable for non-hg19 references.
    merge_bins : bool, optional
        Merge 100kb bins to 5Mb (default ``True``).
    window_size : int, optional
        Merge window size in bp (default 5,000,000).
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    workers : int, optional
        Worker-process count (default 1).
    verbose : int or bool, optional
        Print progress/timing.

    Returns
    -------
    pandas.DataFrame
        DELFI results with the original author's column names.
    """
    if verbose:
        start_time = time.time()
        stderr.write(
            f"""
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
        \n"""
        )
        stderr.write("Reading genome file...\n")

    contigs = chrom_sizes_to_list(chrom_sizes)

    # Resolve the (deprecated) gc_correct flag against no_gc_correct.
    if gc_correct is None:
        gc_correct = not no_gc_correct
    else:
        warnings.warn(
            "Warning: gc_correct is deprecated and may be removed in future "
            "releases. Use no_gc_correct instead"
        )

    gaps = _resolve_gaps(gap_file)

    if verbose:
        stderr.write("Opening bins file...\n")

    bins = pandas.read_csv(
        bins_file,
        names=["contig", "start", "stop"],
        usecols=[0, 1, 2],
        dtype={"contig": str, "start": np.int32, "stop": np.int32},
        delimiter="\t",
        comment="#",
    )

    if verbose:
        stderr.write(f"{bins.shape[0]} bins read from file.\n")
        stderr.write("Filtering gaps...\n")

    if gaps is not None:
        overlaps_gap = overlaps(
            bins["contig"].to_numpy(),
            bins["start"].to_numpy(),
            bins["stop"].to_numpy(),
            gaps.gaps["contig"],
            gaps.gaps["start"],
            gaps.gaps["stop"],
        )
        gapless_bins = bins.loc[~overlaps_gap]
        if verbose:
            stderr.write(
                f"{bins.shape[0] - gapless_bins.shape[0]} bins removed\n"
            )
    else:
        if verbose:
            stderr.write("No gaps specified, skipping.\n")
        gapless_bins = bins

    if verbose:
        stderr.write("Preparing to generate short and long coverages.\n")

    # Pre-parse shared, read-only inputs once. These are handed to each worker
    # via the Pool initializer rather than pickled into every task.
    blacklist_by_contig = _load_blacklist_indexed(blacklist_file)
    contig_gaps_by_contig = {}
    if gaps is not None:
        for contig, _size in contigs:
            contig_gaps_by_contig[contig] = gaps.get_contig_gaps(contig)

    window_args = []
    for contig, size in contigs:
        for _, start, stop, *_ in (
            gapless_bins.loc[gapless_bins.loc[:, "contig"] == contig].itertuples(
                index=False, name=None
            )
        ):
            window_args.append(
                (
                    contig,
                    start,
                    stop,
                    verbose - 1 if verbose > 1 else 0,
                )
            )

    if verbose:
        stderr.write(f"{len(window_args)} windows created.\n")
        stderr.write("Calculating fragment lengths...\n")

    with Pool(
        workers,
        initializer=_delfi_pool_initializer,
        initargs=(
            input_file,
            reference_file,
            quality_threshold,
            blacklist_by_contig,
            contig_gaps_by_contig,
        ),
    ) as pool:
        windows = pool.starmap(_delfi_single_window, tqdm(window_args), 50)

    if verbose:
        stderr.write("Done.\n")
        stderr.write("Removing remaining accrocentric bins...\n")

    window_df = pandas.DataFrame(
        windows,
        columns=[
            "contig",
            "start",
            "stop",
            "arm",
            "short",
            "long",
            "gc",
            "num_frags",
        ],
    )
    trimmed_windows = window_df.loc[window_df["arm"] != "NOARM", :].copy()

    if verbose:
        stderr.write(f"{trimmed_windows.shape[0]} bins remaining...\n")
        stderr.write("Calculating ratio...\n")

    # short/long ratio, guarding against division by zero.
    trimmed_windows["ratio"] = np.where(
        trimmed_windows["long"] == 0,
        np.nan,
        trimmed_windows["short"] / trimmed_windows["long"],
    )

    # Remove the two hg19 no-coverage windows (by positional index).
    if remove_nocov:
        keep = np.logical_and(
            np.arange(trimmed_windows.shape[0]) != 8779,
            np.arange(trimmed_windows.shape[0]) != 13664,
        )
        corrected_delfi_drop_nocov = trimmed_windows.loc[keep].reset_index()
    else:
        corrected_delfi_drop_nocov = trimmed_windows

    if gc_correct:
        if verbose:
            stderr.write("GC bias correction...\n")
        gc_corrected = delfi_gc_correct(
            corrected_delfi_drop_nocov, 0.75, 8, verbose
        )
    else:
        gc_corrected = corrected_delfi_drop_nocov

    if merge_bins:
        if verbose:
            stderr.write("Merging bins...\n")
        final_bins = delfi_merge_bins(gc_corrected, gc_correct, verbose=verbose)
    else:
        final_bins = gc_corrected

    if verbose:
        stderr.write(f"{final_bins.shape[0]} bins remaining.\n")

    if output_file is not None:
        _write_delfi(final_bins, output_file)

    num_frags = sum(window[7] for window in windows)

    if verbose:
        end_time = time.time()
        stderr.write(f"{num_frags} fragments included.\n")
        stderr.write(f"delfi took {end_time - start_time} s to complete\n")
    return final_bins


def _resolve_gaps(gap_file) -> GenomeGaps | None:
    """Coerce the ``gap_file`` argument into a :class:`GenomeGaps` or ``None``."""
    if gap_file is None:
        return None
    if isinstance(gap_file, str):
        return GenomeGaps(gap_file)
    if isinstance(gap_file, GenomeGaps):
        return gap_file
    raise TypeError(f"{type(gap_file)} is not accepted type for gap_file")


def _write_delfi(final_bins: pandas.DataFrame, output_file: str) -> None:
    """Write DELFI results to BED/TSV/CSV/gz or stdout."""
    output_delfi = final_bins.rename(columns={"contig": "#contig"})
    if output_file.endswith(".bed") or output_file.endswith(".tsv"):
        output_delfi.to_csv(output_file, sep="\t", index=False)
    elif output_file.endswith(".csv"):
        final_bins.to_csv(output_file, sep=",", index=False)
    elif output_file.endswith(".bed.gz"):
        output_delfi.to_csv(output_file, sep="\t", index=False, encoding="gzip")
    elif output_file == "-":
        # Stream tab-separated rows to stdout. (The original joined raw,
        # non-string tuple values here and crashed; rows are stringified now.)
        for window in final_bins.itertuples(index=False, name=None):
            stdout.write("\t".join(str(field) for field in window) + "\n")
    else:
        raise ValueError(
            "Invalid file type! Only .bed, .bed.gz, and .tsv suffixes allowed."
        )


def _delfi_single_window(
    contig: str,
    window_start: int,
    window_stop: int,
    verbose: Union[int, bool] = False,
) -> tuple:
    """Compute short/long counts and GC content for one window.

    Reads the alignment, reference, blacklist, and gap inputs from the
    worker-process globals set by ``_delfi_pool_initializer`` so that no file
    is reopened and no input re-parsed per window.
    """
    contig_gaps = (
        _WORKER_CONTIG_GAPS.get(contig)
        if _WORKER_CONTIG_GAPS is not None
        else None
    )

    # Skip windows in centromeres/telomeres or short arms.
    if contig_gaps is not None:
        if contig_gaps.in_tcmere(window_start, window_stop):
            return (contig, window_start, window_stop, "NOARM", np.nan, np.nan, np.nan, 0)
        arm = contig_gaps.get_arm(window_start, window_stop)
        if arm == "NOARM":
            return (contig, window_start, window_stop, "NOARM", np.nan, np.nan, np.nan, 0)
    else:
        arm = contig

    blacklist_regions = _blacklist_in_window(contig, window_start, window_stop)

    short_lengths = 0
    long_lengths = 0
    num_frags = 0

    # Iterate over each fragment in the window using the shared alignment
    # handle. AlignmentWrapper.fetch applies the same quality / read-pair
    # filtering as frag_generator; we replicate frag_generator's length window
    # (100-220) and its default "midpoint" intersect policy here so output is
    # unchanged.
    for frag in _WORKER_ALIGNMENT.fetch(contig, window_start, window_stop):
        frag_start = frag.start
        frag_stop = frag.stop
        frag_length = frag_stop - frag_start

        if frag_length < 100 or frag_length > 220:
            continue

        midpoint = (frag_start + frag_stop) // 2
        if midpoint < window_start or midpoint >= window_stop:
            continue

        blacklisted = False
        for region in blacklist_regions:
            if (
                (frag_start >= region[0] and frag_start < region[1])
                and (frag_stop >= region[0] and frag_stop < region[1])
            ):
                blacklisted = True
                break

        if contig_gaps is not None and contig_gaps.in_tcmere(frag_start, frag_stop):
            continue

        if not blacklisted:
            if frag_length >= 151:
                long_lengths += 1
            else:
                short_lengths += 1
            num_frags += 1

    # GC content of the window from the shared reference handle (vectorized
    # base counting). ref_bases is upper-cased by ReferenceWrapper.
    if not valid_interval(
        _WORKER_REF.chroms, contig, window_start, window_stop, throw_on_error=False
    ):
        logger.warning(
            f"Invalid interval {contig}:{window_start}-{window_stop} for "
            "reference. Skipping GC calculation."
        )
        ref_bases = ""
    else:
        ref_bases = _WORKER_REF.sequence(contig, window_start, window_stop)

    num_gc = ref_bases.count("G") + ref_bases.count("C")

    window_coverage = window_stop - window_start
    gc_content = num_gc / window_coverage if num_frags > 0 else np.nan

    coverage_short = short_lengths
    coverage_long = long_lengths

    if verbose:
        stderr.write(
            f"{contig}:{window_start}-{window_stop} short: "
            f"{coverage_short} long: {coverage_long}, gc_content: "
            f"{gc_content * 100}%\n"
        )

    return (
        contig,
        window_start,
        window_stop,
        arm,
        coverage_short,
        coverage_long,
        gc_content,
        num_frags,
    )
