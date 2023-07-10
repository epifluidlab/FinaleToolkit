from __future__ import annotations

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

def _delfi_gc_adjust(
        windows:np.ndarray,
        verbose:bool=False
):
    """
    Helper function that takes window data and performs GC adjustment.
    """
    #LOESS/LOWESS regression for short and long
    short_loess = delfi_loess(windows['gc'], windows['short'])
    long_loess = delfi_loess(windows['gc'], windows['long'])

    corrected_windows = windows.copy()

    # GC correction
    corrected_windows['short'] = (
        windows['short']
        - short_loess
        + np.nanmedian(windows['short'])
    )

    corrected_windows['long'] = (
        windows['long']
        - long_loess
        + np.nanmedian(windows['long'])
    )

    return corrected_windows


def delfi_loess(gc, coverage):
    """
    Function to simplify loess.
    """
    coverage_loess = lowess(
        coverage,
        gc,
        0.75,
        5,
        return_sorted=False,
        missing='drop'
    )
    return coverage_loess