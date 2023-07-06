from __future__ import annotations
import time
from typing import Union, Tuple
from sys import stdout, stderr
from shutil import get_terminal_size

import numpy as np
import pysam
from tqdm import tqdm
from numba import jit

from finaletools.utils import not_read1_or_low_quality, cli_hist


def single_frag_length(
        input_file: Union[str, pysam.AlignmentFile],
        contig: str=None,
        start: int=None,
        stop: int=None,
        output_file: str=None,
        quality_threshold: int=30,
        verbose: bool=False
    ) -> np.ndarray:
    """
    Return `np.ndarray` containing lengths of fragments in `input_file`
    that are above the quality threshold and are proper-paired reads.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, SAM, or CRAM file containing paired-end fragment reads or
        its path. `AlignmentFile` must be opened in read mode.
    contig : string, optional
    output_file : string, optional
    quality_threshold : int, optional
    verbose : bool, optional

    Returns
    -------
    lengths : numpy.ndarray
        `ndarray` of fragment lengths from file and contig if
        specified.
    """
    if (verbose):
        start_time = time.time()
        stderr.write("Finding frag lengths.\n")

    lengths = []    # list of fragment lengths

    input_is_file = False
    try:
        # handling input types
        if (type(input_file) == pysam.AlignmentFile):
            sam_file = input_file
        elif input_file.endswith('bam'):
            input_is_file = True
            if (verbose):
                stderr.write(f'Opening {input_file}\n')
            sam_file = pysam.AlignmentFile(input_file)
        else:
            raise ValueError(
                'Invalid input_file type. Only BAM or SAM files are allowed.'
            )

        if (verbose):
            stderr.write('Counting reads')
        count = sam_file.count(contig=contig,
                            start=start,
                            stop=stop) if verbose else None
        if (verbose):
            stderr.write(f'{count} reads counted\n')
        # Iterating on each read in file in specified contig/chromosome
        for read1 in (tqdm(sam_file.fetch(contig=contig,
                                        start=start,
                                        stop=stop), total=count)
                    if verbose
                    else sam_file.fetch(contig=contig,
                                        start=start,
                                        stop=stop)):
            # Only select forward strand and filter out non-paired-end
            # reads and low-quality reads
            if (not_read1_or_low_quality(read1, quality_threshold)):
                pass
            else:
                # append length of fragment to list
                lengths.append(abs(read1.template_length))
    finally:
        if input_is_file:
            sam_file.close()

    # convert to array
    lengths = np.array(lengths)

    # check if output specified
    if (type(output_file) == str):
        if output_file.endswith(".bin"): # binary file
            with open(output_file, 'wt') as out:
                lengths.tofile(out)
        elif output_file == '-':
            for line in lengths:
                stdout.write(f'{line}\n')

        else:   # unaccepted file type
            raise ValueError(
                'output_file can only have suffixes .wig or .wig.gz.'
                )

    elif (output_file is not None):
        raise TypeError(
            f'output_file is unsupported type "{type(input_file)}". '
            'output_file should be a string specifying the path of the file '
            'to write output scores to.'
            )

    if (verbose):
        end_time = time.time()
        stderr.write(
            f'frag_length took {end_time - start_time} s to complete\n'
        )

    return lengths


def single_frag_length_bins(
    input_file: Union[str, pysam.AlignmentFile],
    contig: str=None,
    start: int=None,
    stop: int=None,
    bin_size: int=None,
    output_file: str=None,
    contig_by_contig: bool=False,
    histogram: bool=False,
    quality_threshold: int=30,
    verbose: Union[bool, int]=False
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Takes input_file, computes frag lengths of fragments and returns
    two arrays containing bins and counts by size. Options
    for printing data to stdout or file also exist.

    Parameters
    ----------
    input_file : str or AlignmentFile
    contig : str, optional
    start : int, optional
    stop : int, optional
    bin_size : int, optional
    output_file : str, optional
    workers : int, optional

    Returns
    -------
    bins : ndarray
    counts : ndarray
    """
    # generating fragment lengths
    frag_lengths = single_frag_length(
        input_file,
        contig,
        start,
        stop,
        quality_threshold=quality_threshold,
        verbose=verbose-1 if verbose>1 else 0
    )

    # generating bins and counts
    if bin_size is None:
        if histogram:
            term_width, term_height = get_terminal_size((80, 24))
            n_bins = term_width - 12

            start = np.min(frag_lengths)
            stop = np.max(frag_lengths)

            bin_size = round((stop - start) / n_bins)
        else:
            bin_size = 5

    start = np.min(frag_lengths)
    stop = np.max(frag_lengths)
    n_bins = (stop - start) // bin_size

    bins = np.arange(start, stop, bin_size)
    counts = []

    # generate histogram
    for bin in bins:
        count = np.sum((frag_lengths >= bin) * (frag_lengths < (bin + bin_size)))
        counts.append(count)
    bins = np.append(bins, stop)

    if histogram:
        cli_hist(bins, counts)

    # TODO: output to file

    return bins, counts


def interval_frag_length_summary(

):
    return None