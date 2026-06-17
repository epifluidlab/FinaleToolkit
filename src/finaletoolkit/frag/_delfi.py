"""
DELFI fragmentation feature (Cristiano et al., 2019): short/long fragment
ratios per genomic bin, with optional GC correction and 100kb->5Mb merging.
"""
from __future__ import annotations

import time
import warnings
from multiprocessing.pool import Pool
from pathlib import Path
from sys import stderr, stdout
from typing import Union

import numpy as np
import pandas
from tqdm import tqdm

from finaletoolkit.frag._delfi_gc_correct import delfi_gc_correct
from finaletoolkit.frag._delfi_merge_bins import delfi_merge_bins
from finaletoolkit.genome.gaps import ContigGaps, GenomeGaps
from finaletoolkit.io.reference import ReferenceWrapper
from finaletoolkit.utils import chrom_sizes_to_list, frag_generator, overlaps
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

    window_args = []
    for contig, size in contigs:
        contig_gaps = gaps.get_contig_gaps(contig) if gaps is not None else None
        for _, start, stop, *_ in (
            gapless_bins.loc[gapless_bins.loc[:, "contig"] == contig].itertuples(
                index=False, name=None
            )
        ):
            window_args.append(
                (
                    input_file,
                    reference_file,
                    contig_gaps,
                    contig,
                    start,
                    stop,
                    blacklist_file,
                    quality_threshold,
                    verbose - 1 if verbose > 1 else 0,
                )
            )

    if verbose:
        stderr.write(f"{len(window_args)} windows created.\n")
        stderr.write("Calculating fragment lengths...\n")

    with Pool(workers) as pool:
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
    input_file: str,
    reference_file: str | Path,
    contig_gaps: ContigGaps,
    contig: str,
    window_start: int,
    window_stop: int,
    blacklist_file: str = None,
    quality_threshold: int = 30,
    verbose: Union[int, bool] = False,
) -> tuple:
    """Compute short/long counts and GC content for one window."""
    blacklist_regions = []
    if blacklist_file is not None:
        with open(blacklist_file) as blacklist:
            for line in blacklist:
                region_contig, region_start, region_stop, *_ = line.split()
                region_start = int(region_start)
                region_stop = int(region_stop)
                if (
                    contig == region_contig
                    and window_start <= region_start
                    and window_stop >= region_stop
                ):
                    blacklist_regions.append((region_start, region_stop))

    # Skip windows in centromeres/telomeres or short arms.
    if contig_gaps is not None:
        if contig_gaps.in_tcmere(window_start, window_stop):
            return (contig, window_start, window_stop, "NOARM", np.nan, np.nan, np.nan, 0)
        arm = contig_gaps.get_arm(window_start, window_stop)
        if arm == "NOARM":
            return (contig, window_start, window_stop, "NOARM", np.nan, np.nan, np.nan, 0)
    else:
        arm = contig

    short_lengths = 0
    long_lengths = 0
    num_frags = 0

    for _, frag_start, frag_stop, _, _ in frag_generator(
        input_file,
        contig,
        quality_threshold,
        window_start,
        window_stop,
        min_length=100,
        max_length=220,
        reference_file=reference_file,
    ):
        frag_length = frag_stop - frag_start
        assert frag_length > 0, (
            f"Frag length of {frag_length} found at"
            f"{contig}:{frag_start}-{frag_stop}."
        )

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

    # GC content of the window from the reference (vectorized base counting).
    with ReferenceWrapper(reference_file, use_lock=False) as ref_seq:
        if not valid_interval(
            ref_seq.chroms, contig, window_start, window_stop, throw_on_error=False
        ):
            logger.warning(
                f"Invalid interval {contig}:{window_start}-{window_stop} for "
                "reference. Skipping GC calculation."
            )
            ref_bases = ""
        else:
            ref_bases = ref_seq.sequence(contig, window_start, window_stop)

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
