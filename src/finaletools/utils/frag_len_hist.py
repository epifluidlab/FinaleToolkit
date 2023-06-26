from __future__ import annotations
from shutil import get_terminal_size

import numpy as np

def frag_len_hist(
        frag_lengths: np.ndarray
):
    term_width, term_height = get_terminal_size((80, 24))
    n_bins = term_width - 12
    counts, bins = np.histogram(frag_lengths, bins=n_bins)
    max_count = np.max(counts)
    max_height = term_height - 2
    block_height = max_height * 8    # height for utf-8 block elements
    block_counts = np.round(counts / max_count * block_height)

    hist_array = np.zeros((max_height, n_bins), dtype=np.byte)
    for i, block_count in range(n_bins), block_counts:
        solid_height = block_count % 8
        hist_array[-solid_height:] = 8




