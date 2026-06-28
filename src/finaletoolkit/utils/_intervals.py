"""
Helpers for collapsing overlapping intervals read from BED files.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

__all__ = [
    "_merge_overlapping_intervals",
    "_reduce_overlaps_in_file",
    "_convert_to_list",
    "_merge_all_intervals",
]


def _merge_overlapping_intervals(
    intervals: List[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    """Merge overlapping ``(start, stop)`` intervals on a single contig.

    Parameters
    ----------
    intervals : list of (int, int)
        Intervals on one contig. Mutated in place by sorting.

    Returns
    -------
    list of (int, int)
        Non-overlapping intervals sorted by start.
    """
    intervals.sort(key=lambda x: x[0])
    merged: List[Tuple[int, int]] = []
    for interval in intervals:
        if not merged or interval[0] > merged[-1][1]:
            merged.append(interval)
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))
    return merged


def _reduce_overlaps_in_file(interval_file: str) -> Dict[str, List[Tuple[int, int]]]:
    """Read a BED file and merge overlapping intervals per contig.

    Returns
    -------
    dict
        Maps each contig to its list of merged ``(start, stop)`` intervals.
    """
    intervals_dict: Dict[str, List[Tuple[int, int]]] = {}
    with open(interval_file, "r") as file:
        for line in file:
            chrom, start, end = line.strip().split("\t")[:3]
            start, end = int(start), int(end)
            intervals_dict.setdefault(chrom, []).append((start, end))

    return {
        chrom: _merge_overlapping_intervals(intervals)
        for chrom, intervals in intervals_dict.items()
    }


def _convert_to_list(
    reduced_intervals: Dict[str, List[Tuple[int, int]]]
) -> Dict[str, List[list]]:
    """Convert ``{chrom: [(start, stop), ...]}`` to ``{chrom: [[chrom, start, stop], ...]}``."""
    return {
        chrom: [[chrom, start, end] for start, end in intervals]
        for chrom, intervals in reduced_intervals.items()
    }


def _merge_all_intervals(converted_intervals: Dict[str, List[list]]) -> List[list]:
    """Flatten a ``{chrom: [[chrom, start, stop], ...]}`` mapping into one list."""
    all_intervals: List[list] = []
    for intervals in converted_intervals.values():
        all_intervals.extend(intervals)
    return all_intervals
