from __future__ import annotations
from shutil import get_terminal_size

import numpy as np

def cli_hist(
        bins: np.ndarray,
        counts: np.ndarray
):
    term_width, term_height = get_terminal_size((80, 24))
    n_bins = bins.shape[0]

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

