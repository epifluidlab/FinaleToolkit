"""
Merge uniform DELFI bins (100kb by default) into larger windows (5Mb by
default), per chromosome arm.
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pandas as pd

from finaletoolkit.utils.logging import get_logger

logger = get_logger(__name__)

__all__ = ["delfi_merge_bins"]


class _BinLayout(NamedTuple):
    """Result of :func:`_resolve_bin_size`."""

    bin_size: int
    inclusive: bool
    keep: np.ndarray  # boolean mask over the input rows


def _resolve_bin_size(bins: pd.DataFrame, window_size: int) -> _BinLayout:
    """Infer and validate a uniform bin size that evenly divides ``window_size``.

    Both BED-style half-open bins (``stop - start == size``) and end-inclusive
    bins (``stop - start == size - 1`` with each bin starting at the previous
    stop + 1, as in the Cristiano et al. HiC compartment bins) are accepted.
    The last bin of a contig may be shorter than the rest (as produced by e.g.
    ``bedtools makewindows``); such bins are flagged for removal in ``keep``.

    Raises
    ------
    ValueError
        If ``window_size`` is not positive, the bins are empty or not uniform,
        the coordinate convention is mixed, or the bin size does not evenly
        divide ``window_size``.
    """
    if window_size <= 0:
        raise ValueError(f"window_size must be positive, got {window_size}.")
    if bins.shape[0] == 0:
        raise ValueError("Cannot determine bin size: no bins were provided.")

    contigs = bins["contig"].to_numpy()
    starts = bins["start"].to_numpy(dtype=np.int64)
    stops = bins["stop"].to_numpy(dtype=np.int64)
    widths = stops - starts

    # The bin with the greatest start on each contig is its last bin.
    is_last = np.zeros(bins.shape[0], dtype=bool)
    is_last[pd.Series(starts).groupby(contigs).idxmax().to_numpy()] = True

    def _describe(i: int) -> str:
        return f"{contigs[i]}:{starts[i]}-{stops[i]} (width {widths[i]})"

    # Every bin but the last of each contig must share one width. If every
    # contig has a single bin, there is nothing to compare against but each
    # other, so all bins are held to the same width.
    reference = widths[~is_last] if (~is_last).any() else widths
    width = int(reference[0])
    if width <= 0:
        raise ValueError(f"Bin {_describe(0)} has non-positive width.")
    mismatched = np.flatnonzero(
        (widths != width) & ~(is_last & (widths < width))
    )
    if mismatched.size > 0:
        examples = ", ".join(_describe(i) for i in mismatched[:3])
        raise ValueError(
            "Bins must all be the same size to merge them, but found "
            f"{mismatched.size} bin(s) whose width differs from {width}: "
            f"{examples}. Only the last bin of a contig may be shorter. Use "
            "uniform bins or pass merge_bins=False (--no-merge-bins)."
        )
    keep = widths == width

    # Infer the coordinate convention from the gap between adjacent bins.
    order = np.lexsort((starts, contigs))
    same_contig = contigs[order][1:] == contigs[order][:-1]
    gaps = (starts[order][1:] - stops[order][:-1])[same_contig]
    half_open = bool((gaps == 0).any())
    inclusive = bool((gaps == 1).any())
    if half_open and inclusive:
        raise ValueError(
            "Bins mix half-open (next start == previous stop) and "
            "end-inclusive (next start == previous stop + 1) coordinates."
        )
    bin_size = width + 1 if inclusive else width
    logger.info(
        f"Detected {bin_size} bp "
        f"{'end-inclusive' if inclusive else 'half-open'} bins "
        f"(stop - start = {width})."
    )

    dropped = np.flatnonzero(~keep)
    if dropped.size > 0:
        examples = ", ".join(_describe(i) for i in dropped[:3])
        logger.warning(
            f"Dropping {dropped.size} short bin(s) at contig ends: {examples}"
            f"{', ...' if dropped.size > 3 else ''}."
        )

    if window_size < bin_size or window_size % bin_size != 0:
        raise ValueError(
            f"Merge window size ({window_size} bp) must be a positive "
            f"multiple of the bin size ({bin_size} bp). Choose a window size "
            f"such as {max(1, window_size // bin_size) * bin_size} bp."
        )
    return _BinLayout(bin_size, inclusive, keep)


def _aggregate_chunk(chunk: pd.DataFrame, arm: str, include_corrected: bool) -> tuple:
    """Aggregate one chunk of bins into a single merged-bin record."""
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
    *,
    window_size: int = 5000000,
    bin_size: int | None = None,
) -> pd.DataFrame:
    """Merge uniform DELFI bins into non-overlapping windows per arm.

    Each window holds ``window_size // bin_size`` consecutive bins (50 x 100kb
    = 5Mb by default). As in the original DELFI scripts, bins are grouped by
    count, so a window may span more than ``window_size`` bp where bins were
    removed (gaps, blacklist). p-arms are aggregated 5'->3'; q-arms are
    aggregated 3'->5' and then reversed. Partial chunks at arm ends are
    dropped.

    Parameters
    ----------
    hundred_kb_bins : pandas.DataFrame
        Uniform bins with an ``arm`` column and DELFI feature columns. When
        ``gc_corrected`` is ``True`` the ``*_corrected`` columns must be
        present.
    gc_corrected : bool, optional
        Whether the input carries GC-corrected columns to aggregate (default
        ``True``).  Unlike the original implementation this flag is honored, so
        merging works with GC correction disabled.
    verbose : bool, optional
        Unused; kept for signature compatibility.
    window_size : int, optional
        Merged window size in bp (default 5,000,000). Must be a multiple of
        the bin size.
    bin_size : int, optional
        Size of the input bins in bp. If omitted, it is inferred from the
        ``start``/``stop`` columns, accepting half-open or end-inclusive
        coordinates; short bins at contig ends are dropped with a warning.

    Returns
    -------
    pandas.DataFrame
        The merged bins, with the same columns as the input (minus any
        ``index`` column).

    Raises
    ------
    ValueError
        If the bins are not uniform or the bin size does not evenly divide
        ``window_size``.
    """
    if bin_size is None:
        layout = _resolve_bin_size(hundred_kb_bins, window_size)
        hundred_kb_bins = hundred_kb_bins.loc[layout.keep]
        bin_size = layout.bin_size
    elif bin_size <= 0 or window_size < bin_size or window_size % bin_size != 0:
        raise ValueError(
            f"Merge window size ({window_size} bp) must be a positive "
            f"multiple of the bin size ({bin_size} bp)."
        )
    bins_per_window = window_size // bin_size

    merged_bins: list[tuple] = []

    for arm in hundred_kb_bins["arm"].unique():
        arm_bins = hundred_kb_bins[hundred_kb_bins["arm"] == arm].reset_index()
        if "p" in arm:
            for i in range(0, arm_bins.shape[0], bins_per_window):
                chunk = arm_bins.loc[i : i + bins_per_window - 1, :]
                if chunk.shape[0] < bins_per_window:
                    continue
                merged_bins.append(_aggregate_chunk(chunk, arm, gc_corrected))
        elif "q" in arm:
            reversed_bins: list[tuple] = []
            # Stop at -1 (not 0) so the chunk ending at row 0 is reachable
            # when bins_per_window == 1.
            for i in range(arm_bins.shape[0] - 1, -1, -bins_per_window):
                chunk = arm_bins.loc[i - (bins_per_window - 1) : i, :]
                if chunk.shape[0] < bins_per_window:
                    continue
                reversed_bins.append(_aggregate_chunk(chunk, arm, gc_corrected))
            merged_bins.extend(reversed(reversed_bins))
        del arm_bins

    return pd.DataFrame(
        merged_bins,
        columns=hundred_kb_bins.columns[hundred_kb_bins.columns != "index"],
    )
