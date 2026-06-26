"""
GC-bias correction for DELFI windows via LOESS smoothing (Cleveland, 1979).
"""
from __future__ import annotations

import time
from sys import stderr

import numpy as np
import pandas
from loess.loess_1d import loess_1d

__all__ = ["delfi_gc_correct"]

# Columns corrected against GC content.
_GC_CORRECT_COLUMNS = ["short", "long", "num_frags", "ratio"]


def delfi_gc_correct(
    windows: pandas.DataFrame,
    alpha: float = 0.75,
    it: int = 8,
    verbose: bool = False,
) -> pandas.DataFrame:
    """GC-correct DELFI window features with a LOESS fit per feature.

    For each of ``short``, ``long``, ``num_frags`` and ``ratio``, a degree-2
    LOESS curve is fit over GC content and subtracted (re-centered on the
    feature median), producing ``{feature}_corrected`` columns.

    Parameters
    ----------
    windows : pandas.DataFrame
        DELFI windows; must contain a ``gc`` column and the feature columns.
    alpha : float, optional
        LOESS smoothing fraction (default 0.75).
    it : int, optional
        Retained for signature compatibility (unused).
    verbose : bool, optional
        Print timing.

    Returns
    -------
    pandas.DataFrame
        A copy of ``windows`` with ``{feature}_corrected`` columns added.
    """
    if verbose:
        start_time = time.time()
        stderr.write(
            f"""
            windows: {windows}
            alpha: {alpha}
            it: {it}
            verbose: {verbose}
        \n"""
        )

    corrected_windows = windows.copy()
    corrected_windows.replace([np.inf, -np.inf], np.nan, inplace=True)

    valid = corrected_windows.dropna()
    gc_range = np.arange(
        valid["gc"].min(),
        valid["gc"].max() + 0.01,
        0.01,
    )

    loess_lines: dict[str, np.ndarray] = {}
    medians: dict[str, float] = {}
    for column in _GC_CORRECT_COLUMNS:
        _, loess_lines[column], _ = loess_1d(
            valid["gc"].to_numpy(),
            valid[column].to_numpy(),
            xnew=gc_range,
            degree=2,
            frac=alpha,
        )
        medians[column] = valid[column].median()

    for column in _GC_CORRECT_COLUMNS:
        corrected = (
            corrected_windows[column]
            - np.interp(corrected_windows["gc"], gc_range, loess_lines[column])
            + medians[column]
        )
        corrected_windows[f"{column}_corrected"] = corrected

    if verbose:
        end_time = time.time()
        stderr.write(
            f"delfi_gc_correct took {end_time - start_time} s to complete\n"
        )

    return corrected_windows
