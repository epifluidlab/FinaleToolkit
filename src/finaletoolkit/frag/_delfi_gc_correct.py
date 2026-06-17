"""
GC-bias correction for DELFI windows via LOESS smoothing (Cleveland, 1979).
"""
from __future__ import annotations

import time
from sys import stderr, stdin, stdout
from typing import Union

import numpy as np
import pandas
from loess.loess_1d import loess_1d

from finaletoolkit.utils._deprecation import deprecated

__all__ = ["delfi_gc_correct", "cli_delfi_gc_correct"]

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


@deprecated
def cli_delfi_gc_correct(
    input_file: str,
    output_file: str,
    header_lines: int = 1,
    verbose: Union[bool, int] = False,
) -> None:
    """CLI wrapper: GC-correct a raw DELFI BED/TSV and write the result.

    Deprecated: prefer running :func:`delfi` with GC correction enabled.
    """

    def read_stdin():
        readline = stdin.readline()
        while readline:
            yield readline
            readline = stdin.readline()

    if input_file == "-":
        input_file = read_stdin()

    raw_delfi = pandas.read_csv(
        input_file,
        skiprows=1,
        names=[
            "contig",
            "start",
            "stop",
            "arm",
            "short",
            "long",
            "gc",
            "num_frags",
            "ratio",
        ],
        dtype={
            "contig": str,
            "start": np.int32,
            "stop": np.int32,
            "arm": str,
            "short": np.double,
            "long": np.double,
            "gc": np.double,
            "num_frags": np.double,
            "ratio": np.double,
        },
        delimiter="\t",
    )
    corrected_delfi = delfi_gc_correct(raw_delfi)
    output_delfi = corrected_delfi.rename(columns={"contig": "#contig"})

    if verbose:
        stderr.write(
            f"{len(raw_delfi) - output_delfi.shape[0]} bins removed.\n"
        )

    if output_file.endswith(".tsv"):
        output_delfi.to_csv(output_file, sep="\t", index=False)
    elif output_file.endswith(".bed"):
        output_delfi.to_csv(output_file, sep="\t", index=False)
    elif output_file.endswith(".bed.gz"):
        output_delfi.to_csv(output_file, sep="\t", index=False, encoding="gzip")
    elif output_file == "-":
        for window in output_delfi.itertuples(index=False, name=None):
            stdout.write("\t".join(str(field) for field in window) + "\n")
