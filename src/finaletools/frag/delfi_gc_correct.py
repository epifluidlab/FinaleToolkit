from __future__ import annotations
from typing import TextIO, Union
from sys import stdout, stdin
import gzip

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

def delfi_gc_correct(
        windows:np.ndarray,
        verbose:bool=False
):
    """
    Helper function that takes window data and performs GC adjustment.
    """
    #LOESS/LOWESS regression for short and long
    short_loess = _delfi_loess(windows['gc'], windows['short'])
    long_loess = _delfi_loess(windows['gc'], windows['long'])

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


def _delfi_loess(gc, coverage):
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

    raw_delfi = np.genfromtxt(
        input_file,
        dtype=[
            ('contig', '<U32'),
            ('start', 'u8'),
            ('stop', 'u8'),
            ('short', 'f8'),
            ('long', 'f8'),
            ('gc', 'f8')
        ],
        skip_header=header_lines,
    )
    corrected_delfi = delfi_gc_correct(raw_delfi)

    # output
    def _write_out(out: TextIO):
        out.write('#contig\tstart\tstop\tshort\tlong\tgc%\n')
        for window in corrected_delfi:
            out.write(
                f'{window[0]}\t{window[1]}\t{window[2]}\t{window[3]}\t'
                f'{window[4]}\t{window[5]}\t{window[6]}\n')

    if output_file.endswith('.tsv'):
        with open(output_file, 'w') as out:
            out.write('contig\tstart\tstop\tshort\tlong\tgc%\n')
            for window in corrected_delfi:
                out.write(
                    f'{window[0]}\t{window[1]}\t{window[2]}\t{window[3]}\t'
                    f'{window[4]}\t{window[5]}\t{window[6]}\n')
    elif output_file.endswith('.bed'):
        with open(output_file, 'w') as out:
            _write_out(out)
    elif output_file.endswith('.bed.gz'):
        with gzip.open(output_file, 'w') as out:
            _write_out(out)
    elif output_file == '-':
        with stdout as out:
            _write_out(out)
    else:
        raise ValueError(
            'Invalid file type! Only .bed, .bed.gz, and .tsv suffixes allowed.'
        )