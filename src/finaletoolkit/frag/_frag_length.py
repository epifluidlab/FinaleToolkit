from __future__ import annotations
import time
from typing import Union
from sys import stdout, stderr
from multiprocessing import Pool
import gzip
from functools import partial

import numpy as np
import pysam
from tqdm import trange, tqdm
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from finaletoolkit.utils.utils import (
    _get_intervals, frag_generator
)


def plot_histogram(
        data_dict,
        num_bins,
        histogram_path="./frag_length_bins_histogram.png",
        stats=None):
    """
    Generates a histogram from `frag_length_bins` results.
    """
    keys = list(data_dict.keys())
    values = list(data_dict.values())

    fig_size = (6, 4)
    font_size = 12
    plt.figure(figsize=fig_size, dpi=1000)
    plt.hist(keys, bins=num_bins, weights=values, color='salmon',
             edgecolor='white', linewidth=0.1)
    plt.xlabel("Fragment Size (bp)", fontsize=font_size*0.8)
    plt.ylabel("Number of Fragments", fontsize=font_size*0.8)
    plt.xticks(fontsize=font_size * 0.7)
    plt.yticks(fontsize=font_size * 0.7)

    def format_ticks(value, pos):
        if value >= 1e6:
            return '{:1.0f}M'.format(value * 1e-6)
        elif value >= 1e3:
            return '{:1.0f}K'.format(value * 1e-3)
        else:
            return '{:1.0f}'.format(value)

    plt.gca().yaxis.set_major_formatter(FuncFormatter(format_ticks))
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    if stats:
        stats_str = "\n".join([f"{stat[0]}: {stat[1]}" for stat in stats])
        plt.text(
            0.95, 0.95, stats_str, transform=plt.gca().transAxes,
            fontsize=font_size * 0.6, verticalalignment='top',
            horizontalalignment='right',
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    plt.tight_layout()
    plt.savefig(histogram_path)


def _distribution_from_gen(generator):
    value_counts = {}
    for fragment in generator:
        length_of_fragment = fragment[2] - fragment[1]
        if length_of_fragment in value_counts:
            value_counts[length_of_fragment] += 1
        else:
            value_counts[length_of_fragment] = 1
    return value_counts


def _find_median(val_freq_dict):
    val = np.array(list(val_freq_dict.keys()))
    freq = np.array(list(val_freq_dict.values()))
    ord = np.argsort(val)
    val = val[ord]
    freq = freq[ord]
    cdf = np.cumsum(freq)
    
    total_count = cdf[-1]
    if total_count % 2 == 1:
        median_index = np.searchsorted(cdf, total_count // 2)
        median_val = val[median_index]
        return float(median_val)
    else:
        median_indices = np.searchsorted(cdf, [total_count // 2, total_count
                                               // 2 + 1])
        median_vals = val[median_indices]
        median_val = np.mean(median_vals)
        return float(median_val)


def _frag_length_stats(
    input_file: Union[str, pysam.AlignmentFile],
    contig: str,
    start: int,
    stop: int,
    name: str,
    min_length: int,
    max_length: int,
    intersect_policy: str,
    quality_threshold: int,
    verbose: Union[bool, int]
):
    frag_gen = frag_generator(input_file, contig, quality_threshold, start,
                              stop, min_length, max_length, intersect_policy,
                              verbose)
    frag_len_dict = _distribution_from_gen(frag_gen)

    if sum(frag_len_dict.values())==0:
        mean, median, stdev, minimum, maximum = 5*[-1]
    else:
        mean = (sum(value * count for value, count in frag_len_dict.items())
                / sum(frag_len_dict.values()))
        median = _find_median(frag_len_dict)
        variance = sum(count * ((value - mean) ** 2) for value, count
                       in frag_len_dict.items()) / sum(frag_len_dict.values())
        stdev = variance ** 0.5
        minimum = min(frag_len_dict.keys())
        maximum = max(frag_len_dict.keys())

    return contig, start, stop, name, mean, median, stdev, minimum, maximum


def _frag_length_stats_star(partial_frag_stat, interval):
    contig, start, stop, name = interval
    return partial_frag_stat(contig=contig, start=start, stop=stop, name=name)


def frag_length(
    input_file: Union[str, pysam.AlignmentFile, pysam.TabixFile],
    contig: str | None = None,
    start: int | None = None,
    stop: int | None = None,
    intersect_policy: str = "midpoint",
    output_file: str | None = None,
    quality_threshold: int = 30,
    verbose: bool = False
) -> np.ndarray:
    """
    Return `np.ndarray` containing lengths of fragments in `input_file`
    that are above the quality threshold and are proper-paired reads.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, CRAM, or fragment file containing paired-end fragment reads or
        its path. `AlignmentFile` must be opened in read mode.
    contig : string, optional
        Contig or chromosome to get fragments from
    start : int, optional
        0-based left-most coordinate of interval
    stop : int, optional
        1-based right-most coordinate of interval
    intersect_policy : str, optional
        Specifies what policy is used to include fragments in the
        given interval. Default is "midpoint". Policies include:
        - midpoint: the average of end coordinates of a fragment lies
        in the interval.
        - any: any part of the fragment is in the interval.
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
    
    frag_gen = frag_generator(
        input_file=input_file,
        contig=contig,
        quality_threshold=quality_threshold,
        start=start,
        stop=stop,
        min_length=0,
        max_length=1000000000,   #TODO: allow to have None
        intersect_policy=intersect_policy,
        verbose=verbose,
    )

    for contig, frag_start, frag_stop, _, _ in frag_gen:
        lengths.append(frag_stop - frag_start)

    if (verbose):
        stderr.write("Converting to array.\n")

    # convert to array
    lengths = np.array(lengths, dtype=np.int32)

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


def frag_length_bins(
    input_file: Union[str, pysam.AlignmentFile],
    contig: str | None = None,
    start: int | None = None,
    stop: int | None = None,
    min_length: int | None = 0,
    max_length: int | None = None,
    bin_size: int = 1,
    output_file: str | None = None,
    intersect_policy: str = "midpoint",
    quality_threshold: int = 30,
    histogram_path: str | None = None,
    verbose: Union[bool, int] = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Takes input_file, computes frag lengths of fragments and returns
    two arrays containing bins and counts by size. Optionally prints
    data to output as a tab delimited table or histogram.

    Parameters
    ----------
    input_file : str or AlignmentFile
    contig : str, optional
    start : int, optional
    stop : int, optional
    bin_size : int, optional
    output_file : str, optional
    intersect_policy : str, optional
        Specifies what policy is used to include fragments in the
        given interval. Default is "midpoint". Policies include:
        - midpoint: the average of end coordinates of a fragment lies
        in the interval.
        - any: any part of the fragment is in the interval.
    workers : int, optional

    Returns
    -------
    bins : ndarray
    counts : ndarray
    """
    if verbose:
        stderr.write(
            f"""
            input_file: {input_file}
            contig: {contig}
            start: {start}
            stop: {stop}
            bin_size: {bin_size}
            output_file: {output_file}
            intersect_policy: {intersect_policy}
            quality_threshold: {quality_threshold}
            histogram_path: {histogram_path}
            verbose: {verbose}
            \n"""
        )
        start_time = time.time()

    if verbose:
        stderr.write("Generating fragment dictionary. \n")

    frag_gen = frag_generator(
        input_file, contig, quality_threshold, start, stop, min_length,
        max_length, intersect_policy, verbose)
    
    frag_len_dict = _distribution_from_gen(frag_gen)

    mean = (sum(value * count for value, count in frag_len_dict.items())
            / sum(frag_len_dict.values()))
    variance = (sum(count * ((value - mean) ** 2)
                    for value, count in frag_len_dict.items())
                / sum(frag_len_dict.values()))
    
    # get statistics
    stats = []
    stats.append(('mean', mean))
    stats.append(('median', _find_median(frag_len_dict)))
    stats.append(('stdev', variance ** 0.5))
    stats.append(('min', min(frag_len_dict.keys())))
    stats.append(('max', max(frag_len_dict.keys())))

    bin_start = min(frag_len_dict.keys())
    bin_stop = max(frag_len_dict.keys())
    n_bins = (bin_stop - bin_start) // bin_size
    bins = np.arange(bin_start, bin_stop+bin_size, bin_size)

    counts = []

    # generate histogram data
    for i in trange(n_bins+1, disable=not verbose, desc="Binning fragments:"):
        bin_lower = bin_start + i * bin_size
        bin_upper = bin_start + (i + 1) * bin_size
        bin_count = sum(count for length, count in frag_len_dict.items()
                        if bin_lower <= length < bin_upper)
        if bin_count is None:
            bin_count = 0 
        counts.append(bin_count)

    # write results to output
    if output_file is not None:
        try:
            out_is_file = False
            if output_file == '-':
                out = stdout
            elif output_file.endswith('.gz'):
                out_is_file = True
                out = gzip.open(output_file, 'w')
            else:
                out_is_file = True
                out = open(output_file, 'w')

            out.write('min\tmax\tcount\n')
            for bin, count in zip(bins, counts):
                out.write(f'{bin}\t{bin+bin_size-1}\t{count}\n')

            if histogram_path!=None:
                plot_histogram(frag_len_dict, num_bins=n_bins,
                               histogram_path=histogram_path, stats=stats)

        finally:
            if out_is_file:
                out.close()

    # generate histogram figure
    elif histogram_path!=None:
        plot_histogram(frag_len_dict, num_bins=n_bins,
                       histogram_path=histogram_path, stats=stats)

    if verbose:
        stop_time = time.time()
        stderr.write(
            f'frag_length_bins took {stop_time-start_time} s to complete.\n'
        )

    return bins, counts


def frag_length_intervals(
    input_file: Union[str, pysam.AlignmentFile],
    interval_file: str,
    output_file: str | None = None,
    min_length: int | None = 0,
    max_length: int | None = None,
    quality_threshold: int = 30,
    intersect_policy: str = "midpoint",
    workers: int = 1,
    verbose: Union[bool, int] = False,
)->list[tuple[str, int, int, str, float, float, int, int]]:
    """
    Takes fragments from BAM file and calculates fragment length
    statistics for each interval in a BED file. If output specified,
    results will be printed into a tab-delimited file.

    Parameters
    ----------
    input_file : str or AlignmentFile
    interval_file : str
    output_file : str, optional
    quality_threshold : int, optional
    workers : int, optional
    verbose : bool or int, optional

    Returns
    -------
    results: list of (contig, start, stop, name, mean, median, stdev, min, max)'
    """
    if verbose:
        stderr.write(
            f"""
            input_file: {input_file}
            interval_file: {interval_file}
            output_file: {output_file}
            min_length: {min_length}
            max_length: {max_length}
            quality_threshold: {quality_threshold}
            workers: {workers}
            verbose: {verbose}
            \n"""
        )
        start_time = time.time()

    if verbose:
        stderr.write('Creating process pool.\n')
    try:
        pool = Pool(processes=workers)
        if verbose:
            stderr.write('Reading intervals.\n')
        intervals = _get_intervals(interval_file)
        
        partial_frag_stat = partial(
            _frag_length_stats, input_file=input_file,min_length=min_length,
            max_length=max_length, intersect_policy=intersect_policy,
            quality_threshold=quality_threshold, verbose=verbose)
        results = pool.map(partial(_frag_length_stats_star, partial_frag_stat),
                           intervals, chunksize=max(len(intervals)//workers,
                                                    1))
        if verbose:
            tqdm.write('Retrieving fragment statistics for file\n')
        output_is_file = False
        if output_file != None:
            if verbose:
                tqdm.write('Writing results to output. \n')
            try:
                if (output_file.endswith('.bed')
                    or output_file.endswith('.bedgraph')
                ):
                    output_is_file = True
                    output = open(output_file, 'w')
                elif output_file.endswith('.bed.gz'):
                    output = gzip.open(output_file, 'w')
                    output_is_file = True
                elif output_file == '-':
                    output = stdout
                else:
                    raise ValueError(
                        'The output file should have .bed or .bed.gz as as '
                        'suffix.'
                    )
                output.write('contig\tstart\tstop\tname\tmean\tmedian\t'
                             'stdev\tmin\tmax\n')   # type: ignore
                output.write(
                    '\n'.join(
                        '\t'.join(
                            str(element) for element in item)
                        for item in results))   # type: ignore
                output.write('\n')  # type: ignore
            finally:
                if output_is_file:
                    output.close()

    finally:
        pool.close()
        
    if verbose:
        stop_time = time.time()
        stderr.write(
            'Calculating fragment length statistics for intervals took '
            f'{stop_time - start_time} s\n'
        )

    return results