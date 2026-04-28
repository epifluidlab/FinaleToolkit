from __future__ import annotations
import gzip
from pathlib import Path

def _merge_overlapping_intervals(intervals):
    intervals.sort(key=lambda x: x[0])
    merged = []
    for interval in intervals:
        if not merged or interval[0] > merged[-1][1]:
            merged.append(interval)
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))
    return merged


def _reduce_overlaps_in_file(interval_file):
    intervals_dict = {}
    with open(interval_file, 'r') as file:
        for line in file:
            chrom, start, end = line.strip().split('\t')[:3]
            start, end = int(start), int(end)
            if chrom not in intervals_dict:
                intervals_dict[chrom] = []
            intervals_dict[chrom].append((start, end))

    reduced_intervals = {}
    for chrom, intervals in intervals_dict.items():
        reduced_intervals[chrom] = _merge_overlapping_intervals(intervals)
    return reduced_intervals    


def _convert_to_list(reduced_intervals):
    converted_intervals = {}
    for chrom, intervals in reduced_intervals.items():
        converted_intervals[chrom] = [[chrom, start, end] for start, end in intervals]
    return converted_intervals


def _merge_all_intervals(converted_intervals):
    all_intervals = []
    for intervals in converted_intervals.values():
        all_intervals.extend(intervals)
    return all_intervals
