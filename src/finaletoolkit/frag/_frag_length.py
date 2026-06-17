"""
Fragment-length features: raw lengths, binned length distributions, and
per-interval length summary statistics.
"""
from __future__ import annotations

import gzip
import time
import warnings
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from sys import stderr, stdout
from typing import NamedTuple, Union

import numpy as np
import pysam
from tqdm import tqdm

from finaletoolkit.utils import frag_generator, get_intervals
from finaletoolkit.utils.typing import FragFile

__all__ = [
    "frag_length",
    "frag_length_bins",
    "frag_length_intervals",
    "FragLengthStats",
    "plot_histogram",
]


class FragLengthStats(NamedTuple):
    """Per-interval fragment-length summary statistics.

    A drop-in replacement for the original 11-tuple: it unpacks and indexes
    identically while also exposing named fields.  Missing-data intervals use
    ``-1`` for every numeric field, as in the original implementation.
    """

    contig: str
    start: int
    stop: int
    name: str
    mean: float
    median: float
    stdev: float
    minimum: int
    maximum: int
    count: int
    frac_short_reads: float


def plot_histogram(
    data_dict,
    num_bins,
    histogram_path: str = "./frag_length_bins_histogram.png",
    stats=None,
) -> None:
    """Render a fragment-length histogram PNG from ``frag_length_bins`` data.

    Parameters
    ----------
    data_dict : dict
        Mapping of fragment length to count.
    num_bins : int
        Number of histogram bins.
    histogram_path : str, optional
        Output PNG path.
    stats : list of (str, value), optional
        Summary statistics to annotate on the plot.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    keys = list(data_dict.keys())
    values = list(data_dict.values())

    fig_size = (6, 4)
    font_size = 12
    plt.figure(figsize=fig_size, dpi=1000)
    plt.hist(
        keys,
        bins=num_bins,
        weights=values,
        color="salmon",
        edgecolor="white",
        linewidth=0.1,
    )
    plt.xlabel("Fragment Size (bp)", fontsize=font_size * 0.8)
    plt.ylabel("Number of Fragments", fontsize=font_size * 0.8)
    plt.xticks(fontsize=font_size * 0.7)
    plt.yticks(fontsize=font_size * 0.7)

    def format_ticks(value, pos):
        if value >= 1e6:
            return "{:1.0f}M".format(value * 1e-6)
        elif value >= 1e3:
            return "{:1.0f}K".format(value * 1e-3)
        return "{:1.0f}".format(value)

    plt.gca().yaxis.set_major_formatter(FuncFormatter(format_ticks))
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)

    if stats:
        stats_str = "\n".join([f"{stat[0]}: {stat[1]}" for stat in stats])
        plt.text(
            0.95,
            0.95,
            stats_str,
            transform=plt.gca().transAxes,
            fontsize=font_size * 0.6,
            verticalalignment="top",
            horizontalalignment="right",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
        )

    plt.tight_layout()
    plt.savefig(histogram_path)


def _distribution_from_gen(generator) -> dict[int, int]:
    """Count fragments by length from a ``frag_generator`` stream."""
    value_counts: dict[int, int] = {}
    for fragment in generator:
        length_of_fragment = fragment[2] - fragment[1]
        value_counts[length_of_fragment] = value_counts.get(length_of_fragment, 0) + 1
    return value_counts


def _find_median(val_freq_dict: dict[int, int]) -> float:
    """Compute the median of a value->frequency distribution."""
    val = np.array(list(val_freq_dict.keys()))
    freq = np.array(list(val_freq_dict.values()))
    order = np.argsort(val)
    val = val[order]
    freq = freq[order]
    cdf = np.cumsum(freq)

    total_count = cdf[-1]
    if total_count % 2 == 1:
        median_index = np.searchsorted(cdf, total_count // 2)
        return float(val[median_index])
    median_indices = np.searchsorted(
        cdf, [total_count // 2, total_count // 2 + 1]
    )
    return float(np.mean(val[median_indices]))


def _frag_length_stats(
    input_file: FragFile,
    contig: str,
    start: int,
    stop: int,
    name: str,
    min_length: int,
    max_length: int,
    short_reads: int,
    intersect_policy: str,
    quality_threshold: int,
    verbose: Union[bool, int],
    reference_file: str | Path | None = None,
) -> FragLengthStats:
    """Compute fragment-length statistics for a single interval."""
    frag_gen = frag_generator(
        input_file,
        contig,
        quality_threshold,
        start,
        stop,
        min_length,
        max_length,
        intersect_policy,
        verbose,
        reference_file=reference_file,
    )
    frag_len_dict = _distribution_from_gen(frag_gen)

    total_count = sum(frag_len_dict.values())

    if total_count == 0:
        return FragLengthStats(contig, start, stop, name, -1, -1, -1, -1, -1, -1, -1)

    mean = (
        sum(value * count for value, count in frag_len_dict.items()) / total_count
    )
    median = _find_median(frag_len_dict)
    variance = (
        sum(count * ((value - mean) ** 2) for value, count in frag_len_dict.items())
        / total_count
    )
    stdev = variance**0.5
    minimum = min(frag_len_dict.keys())
    maximum = max(frag_len_dict.keys())

    n_short_reads = sum(
        count for length, count in frag_len_dict.items() if length <= short_reads
    )
    frac_short_reads = n_short_reads / total_count

    return FragLengthStats(
        contig,
        start,
        stop,
        name,
        mean,
        median,
        stdev,
        minimum,
        maximum,
        total_count,
        frac_short_reads,
    )


def _frag_length_stats_star(partial_frag_stat, interval) -> FragLengthStats:
    contig, start, stop, name = interval
    return partial_frag_stat(contig=contig, start=start, stop=stop, name=name)


def frag_length(
    input_file: Union[str, pysam.AlignmentFile, pysam.TabixFile],
    contig: str | None = None,
    start: int | None = None,
    stop: int | None = None,
    intersect_policy: str = "midpoint",
    output_file: str | None = None,
    quality_threshold: int = 30,
    verbose: bool = False,
    reference_file: str | Path | None = None,
) -> np.ndarray:
    """Return an array of fragment lengths from an alignment/fragment file.

    Parameters
    ----------
    input_file : str, AlignmentFile, or TabixFile
        BAM, CRAM, or tabix-indexed fragment file (or open pysam handle).
    contig : str, optional
        Restrict to this contig.
    start : int, optional
        0-based left-most coordinate of the interval.
    stop : int, optional
        1-based right-most coordinate of the interval.
    intersect_policy : {"midpoint", "any"}, optional
        Region-membership policy (default ``"midpoint"``).
    output_file : str, optional
        ``.bin`` writes a raw binary array; ``"-"`` writes one length per line
        to stdout.
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    verbose : bool, optional
        Print timing information.
    reference_file : str or Path, optional
        Reference genome (required for CRAM).

    Returns
    -------
    numpy.ndarray
        ``int32`` array of fragment lengths.
    """
    if verbose:
        start_time = time.time()
        stderr.write("Finding frag lengths.\n")

    frag_gen = frag_generator(
        input_file=input_file,
        contig=contig,
        quality_threshold=quality_threshold,
        start=start,
        stop=stop,
        min_length=0,
        max_length=1000000000,
        intersect_policy=intersect_policy,
        verbose=verbose,
        reference_file=reference_file,
    )

    lengths = [frag_stop - frag_start for _, frag_start, frag_stop, _, _ in frag_gen]

    if verbose:
        stderr.write("Converting to array.\n")

    lengths = np.array(lengths, dtype=np.int32)

    if isinstance(output_file, str):
        if output_file.endswith(".bin"):
            with open(output_file, "wt") as out:
                lengths.tofile(out)
        elif output_file == "-":
            for line in lengths:
                stdout.write(f"{line}\n")
        else:
            raise ValueError("output_file can only have suffixes .wig or .wig.gz.")
    elif output_file is not None:
        raise TypeError(
            f'output_file is unsupported type "{type(input_file)}". '
            "output_file should be a string specifying the path of the file "
            "to write output scores to."
        )

    if verbose:
        end_time = time.time()
        stderr.write(f"frag_length took {end_time - start_time} s to complete\n")

    return lengths


def frag_length_bins(
    input_file: FragFile,
    contig: str | None = None,
    start: int | None = None,
    stop: int | None = None,
    min_length: int | None = 0,
    max_length: int | None = None,
    bin_size: int = 1,
    output_file: str | None = None,
    intersect_policy: str = "midpoint",
    quality_threshold: int = 30,
    summary_stats: bool = False,
    short_fraction: int | None = None,
    histogram_path: str | None = None,
    verbose: Union[bool, int] = False,
    reference_file: str | Path | None = None,
) -> tuple[np.ndarray, list]:
    """Bin fragment lengths and optionally write a TSV table or histogram.

    Parameters
    ----------
    input_file : str, AlignmentFile, or TabixFile
        BAM/CRAM/fragment input.
    contig : str, optional
        Restrict to this contig (genome-wide if omitted).
    start, stop : int, optional
        Interval bounds (require ``contig``).
    min_length, max_length : int, optional
        Fragment-length filter applied before binning.
    bin_size : int, optional
        Bin width in bp (default 1).
    output_file : str, optional
        TSV/`.gz` path, or ``"-"`` for stdout.
    intersect_policy : {"midpoint", "any"}, optional
        Region-membership policy.
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    summary_stats : bool, optional
        Append summary statistics as ``#``-comment lines to the TSV.
    short_fraction : int, optional
        If set, add a short-fraction statistic (fragments ``<=`` this length).
    histogram_path : str, optional
        If set, also render a histogram PNG here.
    verbose : bool or int, optional
        Print timing/config information.
    reference_file : str or Path, optional
        Reference genome (required for CRAM).

    Returns
    -------
    bins : numpy.ndarray
        Bin lower bounds.
    counts : list of int
        Fragment count per bin (same length as ``bins``).
    """
    if verbose:
        stderr.write(
            f"""
            input_file: {input_file}
            contig: {contig}
            start: {start}
            stop: {stop}
            bin_size: {bin_size}
            output_file: {output_file}
            intersect_policy: {intersect_policy}
            quality_threshold: {quality_threshold}
            summary_stats: {summary_stats}
            short_fraction: {short_fraction}
            histogram_path: {histogram_path}
            verbose: {verbose}
            \n"""
        )
        start_time = time.time()
        stderr.write("Generating fragment dictionary. \n")

    frag_gen = frag_generator(
        input_file,
        contig,
        quality_threshold,
        start,
        stop,
        min_length,
        max_length,
        intersect_policy,
        verbose,
        reference_file=reference_file,
    )

    frag_len_dict = _distribution_from_gen(frag_gen)

    total_count = sum(frag_len_dict.values())
    if total_count == 0:
        warnings.warn(
            "No fragments found in the specified region. Returning empty result.",
            RuntimeWarning,
            stacklevel=2,
        )
        return np.array([]), np.array([])

    mean = (
        sum(value * count for value, count in frag_len_dict.items()) / total_count
    )
    variance = (
        sum(count * ((value - mean) ** 2) for value, count in frag_len_dict.items())
        / total_count
    )

    stats = [
        ("mean", mean),
        ("median", _find_median(frag_len_dict)),
        ("stdev", variance**0.5),
        ("min", min(frag_len_dict.keys())),
        ("max", max(frag_len_dict.keys())),
        ("total count", total_count),
    ]
    if short_fraction is not None:
        short_coverage = sum(
            count
            for length, count in frag_len_dict.items()
            if length <= short_fraction
        )
        stats.append(
            (f"short fraction (s{short_fraction})", short_coverage / total_count)
        )

    bin_start = min(frag_len_dict.keys())
    bin_stop = max(frag_len_dict.keys())
    n_bins = (bin_stop - bin_start) // bin_size
    bins = np.arange(bin_start, bin_stop + bin_size, bin_size)

    # Vectorized binning: accumulate each length's frequency into its bin.
    lengths_arr = np.fromiter(frag_len_dict.keys(), dtype=np.int64)
    freqs_arr = np.fromiter(frag_len_dict.values(), dtype=np.int64)
    bin_index = (lengths_arr - bin_start) // bin_size
    counts_arr = np.zeros(n_bins + 1, dtype=np.int64)
    np.add.at(counts_arr, bin_index, freqs_arr)
    counts = counts_arr.tolist()

    if output_file is not None:
        out_is_file = False
        try:
            if output_file == "-":
                out = stdout
            elif output_file.endswith(".gz"):
                out_is_file = True
                out = gzip.open(output_file, "wt")
            else:
                out_is_file = True
                out = open(output_file, "w")

            out.write("min\tmax\tcount\n")
            for bin_val, count in zip(bins, counts):
                out.write(f"{bin_val}\t{bin_val + bin_size - 1}\t{count}\n")

            if summary_stats:
                for name, value in stats:
                    out.write(f"#{name}: {value}\n")
        finally:
            if out_is_file:
                out.close()

    if histogram_path is not None:
        plot_histogram(
            frag_len_dict,
            num_bins=n_bins,
            histogram_path=histogram_path,
            stats=stats,
        )

    if verbose:
        stop_time = time.time()
        stderr.write(
            f"frag_length_bins took {stop_time - start_time} s to complete.\n"
        )

    return bins, counts


def frag_length_intervals(
    input_file: Union[str, pysam.AlignmentFile],
    interval_file: str,
    output_file: str | None = None,
    min_length: int | None = 0,
    max_length: int | None = None,
    quality_threshold: int = 30,
    intersect_policy: str = "midpoint",
    short_reads: int = 150,
    workers: int = 1,
    verbose: Union[bool, int] = False,
    reference_file: str | Path | None = None,
) -> list[FragLengthStats]:
    """Compute per-interval fragment-length statistics over a BED file.

    Parameters
    ----------
    input_file : str or AlignmentFile
        BAM/CRAM/fragment input.
    interval_file : str
        BED file of intervals.
    output_file : str, optional
        BED/`.gz` path or ``"-"`` for stdout.
    min_length, max_length : int, optional
        Fragment-length filter.
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    intersect_policy : {"midpoint", "any"}, optional
        Region-membership policy.
    short_reads : int, optional
        Short-read length cutoff for the short fraction (default 150).
    workers : int, optional
        Worker-process count (default 1).
    verbose : bool or int, optional
        Print timing/config information.
    reference_file : str or Path, optional
        Reference genome (required for CRAM).

    Returns
    -------
    list of FragLengthStats
        One record per interval (``contig, start, stop, name, mean, median,
        stdev, min, max, count, frac_short_reads``).
    """
    if verbose:
        stderr.write(
            f"""
            input_file: {input_file}
            interval_file: {interval_file}
            output_file: {output_file}
            min_length: {min_length}
            max_length: {max_length}
            quality_threshold: {quality_threshold}
            workers: {workers}
            verbose: {verbose}
            \n"""
        )
        start_time = time.time()
        stderr.write("Creating process pool.\n")

    pool = Pool(processes=workers)
    try:
        if verbose:
            stderr.write("Reading intervals.\n")
        intervals = get_intervals(interval_file)

        partial_frag_stat = partial(
            _frag_length_stats,
            input_file=input_file,
            min_length=min_length,
            max_length=max_length,
            short_reads=short_reads,
            intersect_policy=intersect_policy,
            quality_threshold=quality_threshold,
            verbose=verbose,
            reference_file=reference_file,
        )

        results = pool.map(
            partial(_frag_length_stats_star, partial_frag_stat),
            intervals,
            chunksize=max(len(intervals) // workers, 1),
        )

        if verbose:
            tqdm.write("Retrieving fragment statistics for file\n")

        output_is_file = False
        if output_file is not None:
            if verbose:
                tqdm.write("Writing results to output. \n")
            try:
                if output_file.endswith(".bed") or output_file.endswith(".bedgraph"):
                    output_is_file = True
                    output = open(output_file, "w")
                elif output_file.endswith(".bed.gz"):
                    output = gzip.open(output_file, "w")
                    output_is_file = True
                elif output_file == "-":
                    output = stdout
                else:
                    raise ValueError(
                        "The output file should have .bed or .bed.gz as as suffix."
                    )
                output.write(
                    "contig\tstart\tstop\tname\tmean\tmedian\t"
                    "stdev\tmin\tmax\tcount"
                    f"\ts{short_reads}\n"
                )
                output.write(
                    "\n".join(
                        "\t".join(str(element) for element in item)
                        for item in results
                    )
                )
                output.write("\n")
            finally:
                if output_is_file:
                    output.close()
    finally:
        pool.close()

    if verbose:
        stop_time = time.time()
        stderr.write(
            "Calculating fragment length statistics for intervals took "
            f"{stop_time - start_time} s\n"
        )

    return results
