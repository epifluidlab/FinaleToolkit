from __future__ import annotations
from typing import TextIO, Union
from sys import stdout, stdin
import gzip
import time
from sys import stderr

import pandas
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
from loess.loess_1d import loess_1d


def delfi_gc_correct(
        windows: pandas.DataFrame,
        alpha: float = 0.75,
        it: int = 8,
        verbose:bool=False
):
    """
    Helper function that takes window data and performs GC adjustment.
    """
    
    if (verbose):
        start_time = time.time()
        stderr.write(f"""
            windows: {windows}
            alpha: {alpha}
            it: {it}
            verbose: {verbose}
        \n""")

    ##### loess fit with linear interpretation for performance #####
        
    corrected_windows = windows.copy()  # avoid inplace modification
    corrected_windows.replace([np.inf, -np.inf], np.nan, inplace=True)
    gc_range = np.arange(   # range over gc values to save time
        corrected_windows.dropna()["gc"].min(),
        corrected_windows.dropna()["gc"].max()+0.01,
        0.01)
    loess_lines = dict()
    medians = dict()
    column_names = ['short', 'long', 'num_frags', 'ratio']

    # calculate loess curve over range of gc values and store to dict
    for column in column_names:
        _, loess_lines[column], _ = loess_1d(
            corrected_windows.dropna()["gc"].to_numpy(),
            corrected_windows.dropna()[column].to_numpy(),
            xnew=gc_range,
            degree=2,
            frac=0.75)

        medians[column] = corrected_windows.dropna()[column].median()

    # adding new columns
    for column in column_names:
        corrected = (
            corrected_windows[column]
            - np.interp(corrected_windows['gc'],
                        gc_range,
                        loess_lines[column]) + medians[column]
            )
        corrected_windows[f'{column}_corrected'] = corrected

    if (verbose):
        end_time = time.time()
        stderr.write(f'delfi_gc_correct took {end_time - start_time} s to complete\n')

    return corrected_windows


def cli_delfi_gc_correct(
        input_file: str,
        output_file: str,
        header_lines: int=1,
        verbose: Union[bool, int]=False,
):
    """
    Takes path to a BED3+3 containing short and long fractions for DELFI
    and performs GC correction using a LOESS smoother (Cleveland, 1979).
    Smoothing is performed on short and long fractions separately.
    """
    def read_stdin():
        readline = stdin.readline()
        while readline:
            yield readline
            readline = stdin.readline()

    if input_file == '-':
        input_file = read_stdin()

    raw_delfi = pandas.read_csv(
        input_file,
        skiprows=1,
        names=["contig","start","stop","arm","short","long","gc","num_frags","ratio"],
        dtype={"contig":str, "start":np.int32, "stop":np.int32, "arm":str, "short":np.double, "long":np.double, "gc":np.double, "num_frags":np.double, "ratio":np.double,},
        delimiter='\t',
)
    corrected_delfi = delfi_gc_correct(raw_delfi)
    output_delfi = corrected_delfi.rename(columns={'contig':'#contig'})

    # output
    if (verbose):
        stderr.write(f'{len(raw_delfi)-output_delfi.shape[0]} bins '
                     'removed.\n')
        
    if output_file.endswith('.tsv'):
        output_delfi.to_csv(output_file, sep='\t', index=False)
    elif output_file.endswith('.bed'):
        output_delfi.to_csv(
            output_file,
            sep='\t',
            index=False)
    elif output_file.endswith('.bed.gz'):
        output_delfi.to_csv(
            output_file,
            sep='\t',
            index=False,
            encoding='gzip')
    elif output_file == '-':
        with stdout as out:
            for window in output_delfi.itertuples():
                tab_sep = "\t".join(window)
                out.write(
                    f'{tab_sep}\n')