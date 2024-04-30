from __future__ import annotations
from shutil import get_terminal_size
from typing import TextIO
from sys import stdout

import numpy as np

def _cli_hist(
        bins: np.ndarray,
        counts: np.ndarray,
        n_bins: int,
        stats: list,
        out: TextIO,
        title: str=None,
):
    term_width, term_height = get_terminal_size((80, 24))

    max_count = np.max(counts)
    total_count = np.sum(counts)
    max_height = term_height - 5
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

    out.write(f'{title}\n')

    stat_len = len(stats)
    for i in range(max_height):
        row = hist_array[i]
        out.write(f'{((max_height - i) / max_height) * max_count/total_count*100:0>5.2f}%   ')
        out.write(''.join([bars[num] for num in row]))
        if i < stat_len:
            stat = stats[i]
            out.write(f' {(stat[0]+"          ")[:10]}:{stat[1]:02.2f}\n')
        else:
            out.write('\n')
    out.write('len (nt)')
    out.write('   '.join([f'{loc:0>3d}' for loc in bins[0::6]]))
    out.write('\n')
