"""
Post-process raw WPS bigWig signal: local median/mean filtering plus an
optional Savitzky-Golay filter (Savitzky & Golay, 1964).
"""
from __future__ import annotations

import gzip
import traceback
from multiprocessing import Pool
from sys import stderr
from time import time
from typing import Union

import numpy as np
import pyBigWig as pbw
from numpy.lib.stride_tricks import sliding_window_view
from numpy.typing import NDArray
from scipy.signal import savgol_filter

from finaletoolkit.utils import chrom_sizes_to_list

__all__ = ["adjust_wps"]


def _running_stat(data: NDArray, window_size: int, use_mean: bool) -> NDArray:
    """Running median/mean over the first ``len(data) - window_size`` windows.

    Vectorized via :func:`sliding_window_view`; produces values identical to the
    original per-index loop.
    """
    n_windows = len(data) - window_size
    if n_windows <= 0:
        return np.array([], dtype=np.float64)
    windows = sliding_window_view(data, window_size)[:n_windows]
    return np.mean(windows, axis=1) if use_mean else np.median(windows, axis=1)


def _local_filter(
    positions: NDArray, data: NDArray, window_size: int, use_mean: bool
):
    """Subtract a running median/mean from ``data`` and trim ``positions``."""
    running = _running_stat(data, window_size, use_mean)
    adjusted_data = data[window_size // 2 : -(window_size // 2)] - running
    adjusted_positions = positions[window_size // 2 : -(window_size // 2)]
    return adjusted_positions, adjusted_data


def _median_filter(positions: NDArray, data: NDArray, window_size: int):
    """Locally-adjusted running median."""
    return _local_filter(positions, data, window_size, use_mean=False)


def _mean_filter(positions: NDArray, data: NDArray, window_size: int):
    """Locally-adjusted running mean."""
    return _local_filter(positions, data, window_size, use_mean=True)


def _single_adjust_wps_star(args):
    """Unpack a tuple of args and call :func:`_single_adjust_wps`."""
    return _single_adjust_wps(*args)


def _single_adjust_wps(
    input_file: str,
    contig: str,
    start: int,
    stop: int,
    median_window_size: int,
    savgol_window_size: int,
    savgol_poly_deg: int,
    mean: bool,
    subtract_edges: bool,
    edge_size: int,
    savgol: bool,
):
    """Filter the raw WPS of one interval; returns contigs/positions/stops/values."""
    raw_wps = None
    try:
        if input_file.endswith(".bw"):
            raw_wps = pbw.open(input_file, "r")
        else:
            raise ValueError("Invalid filetype for input_file.")

        genomic_range = raw_wps.intervals(contig, start, stop)
        if genomic_range is None:
            # No coverage in this region (or the bigWig was filtered).
            adjusted_positions = np.zeros((0,), dtype=np.int64)
            stops = np.zeros((0,), dtype=np.int64)
            filtered_scores = np.zeros((0,), dtype=np.float64)
            raw_wps.close()
            stderr.write(
                f"No entries in range: {contig}:{start}-{stop}. "
                "This interval will be skipped.\n"
            )
            return (
                len(adjusted_positions) * [contig],
                adjusted_positions,
                stops,
                filtered_scores,
            )

        intervals = np.array(
            list(genomic_range),
            dtype=[("starts", "<i8"), ("stops", "<i8"), ("scores", "<f8")],
        )

        if not all(
            pos1 + 1 == pos2
            for pos1, pos2 in zip(
                intervals["starts"][:-1], intervals["starts"][1:]
            )
        ):
            raise ValueError(
                "BigWig was found to be nonsequential. There may be multiple "
                "entries for one position or gaps in the regions specified in "
                "the interval file."
            )

        if subtract_edges:
            start_mean = np.mean(intervals["scores"][:edge_size])
            stop_mean = np.mean(intervals["scores"][-edge_size:])
            mean_val = np.mean([start_mean, stop_mean])
            intervals["scores"] = intervals["scores"] - mean_val

        if median_window_size > intervals["scores"].shape[0]:
            raise ValueError(
                f"median_window_size ({median_window_size}) cannot be greater "
                f"than the length of interval ({intervals['scores'].shape[0]})."
            )

        adjusted_positions, adjusted_scores = _local_filter(
            intervals["starts"], intervals["scores"], median_window_size, mean
        )

        if savgol:
            filtered_scores = savgol_filter(
                adjusted_scores, savgol_window_size, savgol_poly_deg
            )
        else:
            filtered_scores = adjusted_scores

        assert len(adjusted_positions) == len(filtered_scores)
        stops = adjusted_positions + 1

    except RuntimeError as e:
        traceback.print_exception(e)
        stderr.write(
            "Invalid interval detected:\n"
            f"{contig}:{start}-{stop}. This interval will be skipped.\n"
        )
        adjusted_positions = np.zeros((0,), dtype=np.int64)
        stops = np.zeros((0,), dtype=np.int64)
        filtered_scores = np.zeros((0,), dtype=np.float64)
    finally:
        if raw_wps is not None:
            raw_wps.close()

    return (
        len(adjusted_positions) * [contig],
        adjusted_positions,
        stops,
        filtered_scores,
    )


def adjust_wps(
    input_file: str,
    interval_file: str,
    output_file: str,
    chrom_sizes: str,
    interval_size: int = 5000,
    median_window_size: int = 1000,
    savgol_window_size: int = 21,
    savgol_poly_deg: int = 2,
    savgol: bool = True,
    mean: bool = False,
    subtract_edges: bool = False,
    edge_size: int = 500,
    workers: int = 1,
    verbose: Union[bool, int] = False,
) -> None:
    """Adjust raw WPS in a bigWig with median/mean and Savitzky-Golay filters.

    Parameters
    ----------
    input_file : str
        Path to a bigWig of raw WPS.
    interval_file : str
        BED of intervals WPS was computed over.
    output_file : str
        Output bigWig path.
    chrom_sizes : str
        ``.chrom.sizes`` file for the bigWig header.
    interval_size : int, optional
        Size of each centered interval (default 5000).
    median_window_size : int, optional
        Median/mean filter window (default 1000).
    savgol_window_size : int, optional
        Savitzky-Golay window (default 21).
    savgol_poly_deg : int, optional
        Savitzky-Golay polynomial degree (default 2).
    savgol : bool, optional
        Apply the Savitzky-Golay filter (default ``True``).
    mean : bool, optional
        Use a mean filter instead of median (default ``False``).
    subtract_edges : bool, optional
        Subtract the mean of the interval edges before filtering.
    edge_size : int, optional
        Edge width for ``subtract_edges`` (default 500).
    workers : int, optional
        Worker-process count (default 1).
    verbose : bool or int, optional
        Print progress/timing.
    """
    if verbose:
        start_time = time()
        stderr.write("Reading intervals from bed...\n")

    left_of_site = round(-interval_size / 2)
    right_of_site = round(interval_size / 2)
    assert right_of_site - left_of_site == interval_size

    if not (interval_file.endswith(".bed") or interval_file.endswith(".bed.gz")):
        raise ValueError("Invalid filetype for interval_file.")

    # Amount removed from each end by the median filter; used to detect overlap.
    end_decrease = median_window_size // 2
    intervals = []
    opener = gzip.open if interval_file.endswith(".gz") else open
    with opener(interval_file, "rt") as file:
        for line in file:
            contents = line.split("\t")
            contig = contents[0].strip()
            midpoint = (int(contents[1]) + int(contents[2])) // 2

            start = max(0, midpoint + int(left_of_site))
            stop = midpoint + int(right_of_site)

            # Merge with the previous interval if they would overlap after the
            # median filter shrinks each end (avoids duplicate bigWig entries).
            if (
                len(intervals) > 0
                and intervals[-1][1] == contig
                and intervals[-1][3] - end_decrease > start + end_decrease
            ):
                start = intervals[-1][2]
                intervals.pop(-1)

            intervals.append(
                (
                    input_file,
                    contig,
                    int(start),
                    int(stop),
                    median_window_size,
                    savgol_window_size,
                    savgol_poly_deg,
                    mean,
                    subtract_edges,
                    edge_size,
                    savgol,
                )
            )

    if verbose:
        stderr.write("Opening pool...\n")

    pool = Pool(workers)
    try:
        processed_scores = pool.imap(_single_adjust_wps_star, intervals)

        if verbose:
            stderr.write("Writing to output\n")

        with pbw.open(output_file, "w") as output_bw:
            output_bw.addHeader(chrom_sizes_to_list(chrom_sizes))
            for scores in processed_scores:
                contigs, starts, stops, values = scores
                if len(contigs) == 0:
                    continue
                try:
                    output_bw.addEntries(
                        contigs, starts, ends=stops, values=values
                    )
                except RuntimeError as e:
                    traceback.print_exception(e)
                    stderr.write(
                        "RuntimeError encountered while writing to "
                        f"{output_file} at interval {contigs[0]}:"
                        f"{starts[0]}-{stops[-1]}\n"
                    )
    finally:
        pool.close()

    if verbose:
        end_time = time()
        stderr.write(f"Adjust-WPS took {end_time - start_time} s to run.\n")
