from __future__ import annotations
import gzip
import numpy as np
import pysam
from multiprocessing import Pool
from typing import Union, Tuple
from tqdm import tqdm
from sys import stdout, stderr
import time
from finaletools.utils.utils import _get_intervals, frag_generator
from collections import Counter
import itertools
import matplotlib.pyplot as plt
from functools import partial
from matplotlib.ticker import FuncFormatter

def plot_histogram(data_dict, num_bins, histogram_path="./frag_length_bins_histogram.png", stats=None):
    keys = list(data_dict.keys())
    values = list(data_dict.values())

    fig_size = (6, 4)
    font_size = 12
    plt.figure(figsize=fig_size, dpi=1000)
    plt.hist(keys, bins=num_bins, weights=values, color='salmon', edgecolor='white', linewidth=0.1)
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
        plt.text(0.95, 0.95, stats_str, transform=plt.gca().transAxes,
                 fontsize=font_size * 0.6, verticalalignment='top', horizontalalignment='right',
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
        median_indices = np.searchsorted(cdf, [total_count // 2, total_count // 2 + 1])
        median_vals = val[median_indices]
        median_val = np.mean(median_vals)
        return float(median_val)

    

def frag_length_bins(
    input_file: Union[str, pysam.AlignmentFile],
    contig: str,
    start: int,
    stop: int,
    min_length: int,
    max_length: int,
    bin_size: int,
    intersect_policy: str,
    output_file: str,
    quality_threshold: int,
    histogram_path: str,
    verbose: bool
) -> Tuple[np.ndarray, np.ndarray]:
    if verbose:
        stderr.write(
            f"""
            input_file: {input_file}
            contig: {contig}
            start: {start}
            stop: {stop}
            min_length: {min_length}
            max_length: {max_length}
            bin_size: {bin_size}
            intersect_policy: {intersect_policy}
            output_file: {output_file}
            quality_threshold: {quality_threshold}
            histogram_path: {histogram_path}
            verbose: {verbose}
            \n"""
        )
        start_time = time.time()

    if verbose:
        stderr.write("Generating fragment dictionary. \n")
    frag_gen = frag_generator(input_file, contig, quality_threshold, start, stop, min_length, max_length, intersect_policy, verbose)
    frag_len_dict=_distribution_from_gen(frag_gen)
    mean = sum(value * count for value, count in frag_len_dict.items()) / sum(frag_len_dict.values())
    variance = sum(count * ((value - mean) ** 2) for value, count in frag_len_dict.items()) / sum(frag_len_dict.values())

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

    for i in tqdm(range(n_bins+1), disable=not verbose, desc="Binning fragments..."):
        bin_lower = bin_start + i * bin_size
        bin_upper = bin_start + (i + 1) * bin_size
        bin_count = sum(count for length, count in frag_len_dict.items() if bin_lower <= length < bin_upper)
        counts.append(bin_count)

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
                out.write(f'{bin}\t{bin+bin_size}\t{count}\n')
                
            if histogram_path!=None:
                plot_histogram(frag_len_dict, num_bins=n_bins, histogram_path=histogram_path, stats=stats)
                
        finally:
            if out_is_file:
                out.close()
    elif histogram_path!=None:
        plot_histogram(frag_len_dict, num_bins=n_bins, histogram_path=histogram_path, stats=stats)

    if verbose:
        stop_time = time.time()
        stderr.write(
            f'Calculating fragment length bins took {stop_time-start_time} s to complete.\n'
        )

    return bins, counts



def _frag_length_stats(
    input_file: Union(str, pysam.AlignmentFile),
    contig: str,
    start: int,
    stop: int,
    name: str,
    min_length: int,
    max_length: int,
    intersect_policy: str,
    quality_threshold: int,
    verbose: Union(bool, int)
):
    frag_gen = frag_generator(input_file, contig, quality_threshold, start, stop, min_length, max_length, intersect_policy, verbose)
    frag_len_dict=_distribution_from_gen(frag_gen)

    if sum(frag_len_dict.values())==0:
        mean, median, stdev, minimum, maximum = 5*[-1]
    else:
        mean = sum(value * count for value, count in frag_len_dict.items()) / sum(frag_len_dict.values())
        median = _find_median(frag_len_dict)
        variance = sum(count * ((value - mean) ** 2) for value, count in frag_len_dict.items()) / sum(frag_len_dict.values())
        stdev = variance ** 0.5
        minimum = min(frag_len_dict.keys())
        maximum = max(frag_len_dict.keys())

    return contig, start, stop, name, mean, median, stdev, minimum, maximum


def _frag_length_stats_star(partial_frag_stat, interval):
    contig, start, stop, name = interval
    return partial_frag_stat(contig=contig, start=start, stop=stop, name=name)


def frag_length_intervals(
    input_file: Union[str, pysam.AlignmentFile],
    interval_file: str,
    output_file: str,
    min_length: int,
    max_length: int,
    intersect_policy: str,
    quality_threshold: int,
    workers: int,
    verbose: bool
):
    if verbose:
        stderr.write(
            f"""
            input_file: {input_file}
            interval_file: {interval_file}
            output_file: {output_file}
            min_length: {min_length}
            max_length: {max_length}
            intersect_policy: {intersect_policy}
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
        
        partial_frag_stat = partial(_frag_length_stats, input_file=input_file, min_length=min_length, max_length=max_length, intersect_policy=intersect_policy, quality_threshold=quality_threshold, verbose=verbose)
        iter_results = pool.imap(partial(_frag_length_stats_star, partial_frag_stat), intervals, chunksize=max(len(intervals)//workers, 1))

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
                        'The output file should have .bed or .bed.gz as as suffix.'
                    )
                output.write('contig\tstart\tstop\tname\tmean\tmedian\tstdev\tmin\tmax\n')
                output.write('\n'.join('\t'.join(str(element) for element in item) for item in iter_results))
                output.write('\n')


            finally:
                if output_is_file:
                    output.close()

        results = [result for result in iter_results]

    finally:
        pool.close()

    if verbose:
        stop_time = time.time()
        stderr.write(
            f'Calculating fragment length statistics for intervals took {stop_time-start_time} s\n'
        )

    return #results
