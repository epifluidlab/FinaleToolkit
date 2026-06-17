"""
Merge 100kb DELFI bins into 5Mb (50-bin) windows, per chromosome arm.
"""
from __future__ import annotations

import pandas as pd

__all__ = ["delfi_merge_bins"]

_BINS_PER_WINDOW = 50


def _aggregate_chunk(chunk: pd.DataFrame, arm: str, include_corrected: bool) -> tuple:
    """Aggregate one 50-bin chunk into a single 5Mb-bin record."""
    contig = arm[:-1]
    record = [
        contig,
        chunk["start"].min(),
        chunk["stop"].max(),
        arm,
        chunk["short"].sum(),
        chunk["long"].sum(),
        chunk["gc"].mean(),
        chunk["num_frags"].sum(),
        chunk["ratio"].mean(),
    ]
    if include_corrected:
        record.extend(
            [
                chunk["short_corrected"].sum(),
                chunk["long_corrected"].sum(),
                chunk["num_frags_corrected"].sum(),
                chunk["ratio_corrected"].mean(),
            ]
        )
    return tuple(record)


def delfi_merge_bins(
    hundred_kb_bins: pd.DataFrame,
    gc_corrected: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    """Merge 100kb DELFI bins into non-overlapping 5Mb windows per arm.

    p-arms are aggregated 5'->3'; q-arms are aggregated 3'->5' and then
    reversed, matching the original DELFI scripts.  Partial (<50-bin) chunks at
    arm ends are dropped.

    Parameters
    ----------
    hundred_kb_bins : pandas.DataFrame
        100kb bins with an ``arm`` column and DELFI feature columns.  When
        ``gc_corrected`` is ``True`` the ``*_corrected`` columns must be
        present.
    gc_corrected : bool, optional
        Whether the input carries GC-corrected columns to aggregate (default
        ``True``).  Unlike the original implementation this flag is honored, so
        merging works with GC correction disabled.
    verbose : bool, optional
        Unused; kept for signature compatibility.

    Returns
    -------
    pandas.DataFrame
        The merged 5Mb bins, with the same columns as the input (minus any
        ``index`` column).
    """
    five_mb_bins: list[tuple] = []

    for arm in hundred_kb_bins["arm"].unique():
        arm_bins = hundred_kb_bins[hundred_kb_bins["arm"] == arm].reset_index()
        if "p" in arm:
            for i in range(0, arm_bins.shape[0], _BINS_PER_WINDOW):
                chunk = arm_bins.loc[i : i + _BINS_PER_WINDOW - 1, :]
                if chunk.shape[0] < _BINS_PER_WINDOW:
                    continue
                five_mb_bins.append(_aggregate_chunk(chunk, arm, gc_corrected))
        elif "q" in arm:
            reversed_bins: list[tuple] = []
            for i in range(arm_bins.shape[0] - 1, 0, -_BINS_PER_WINDOW):
                chunk = arm_bins.loc[i - (_BINS_PER_WINDOW - 1) : i, :]
                if chunk.shape[0] < _BINS_PER_WINDOW:
                    continue
                reversed_bins.append(_aggregate_chunk(chunk, arm, gc_corrected))
            five_mb_bins.extend(reversed(reversed_bins))
        del arm_bins

    return pd.DataFrame(
        five_mb_bins,
        columns=hundred_kb_bins.columns[hundred_kb_bins.columns != "index"],
    )
