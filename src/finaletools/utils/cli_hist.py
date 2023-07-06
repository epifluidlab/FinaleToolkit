from __future__ import annotations
from shutil import get_terminal_size

import numpy as np

def frag_len_hist(
        frag_lengths: np.ndarray
):
    term_width, term_height = get_terminal_size((80, 24))
    n_bins = term_width - 12

    start = np.min(frag_lengths)
    stop = np.max(frag_lengths)

    bin_size = round((stop - start) / n_bins)

    print(bin_size)

    n_bins = (stop - start) // bin_size

    bins = np.arange(start, stop, bin_size)
    counts = []
    # generate histogram
    for bin in bins:
        count = np.sum((frag_lengths >= bin) * (frag_lengths < (bin + bin_size)))
        counts.append(count)
    bins = np.append(bins, stop)

    max_count = np.max(counts)
    max_height = term_height - 2
    block_height = max_height * 8 - 1    # height for utf-8 block elements
    block_counts = np.rint(counts / max_count * block_height)

    hist_array = np.zeros((max_height, n_bins), dtype=int)
    for i in range(n_bins):
        block_count = block_counts[i]
        solid_height = int(block_count // 8)
        if solid_height != 0:
            hist_array[-solid_height:, i] = 8
        hist_array[-solid_height-1, i] = block_count % 8

    bars = {
        0: ' ',
        1: '\u2581',
        2: '\u2582',
        3: '\u2583',
        4: '\u2584',
        5: '\u2585',
        6: '\u2586',
        7: '\u2587',
        8: '\u2588'
    }

    for row in hist_array:
        print(''.join([bars[num] for num in row]))

