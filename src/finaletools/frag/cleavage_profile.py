"""Cleavage Profiler

This module is used to find the cleavage profile as described in Zhou
et al 2022 (https://doi.org/10.1073/pnas.2209852119). Cleavage profile
describes the proportion of fragment ends at a site over the depth at
the site (cleavage proportion) calculated over a 5+/- nt window around a
CpG site.
"""

from __future__ import annotations
from typing import Union

import numpy as np

from finaletools.utils.utils import frag_array, overlaps


def cleavage_profile(
    input_file: str,
    contig: str,
    start: int,
    stop: int,
    fraction_low: int=None,
    fraction_high: int=None,
    quality_threshold: int=30,
    verbose: Union[bool, int]=0
) -> np.ndarray:
    frags = frag_array(
        input_file=input_file,
        contig=contig,
        quality_threshold=quality_threshold,
        start=start,
        stop=stop,
        fraction_low=fraction_low,
        fraction_high=fraction_high
    )

    profiles = np.zeros(stop-start, dtype=[
        ('pos', 'i8'),
        ('ends', 'i8'),
        ('overlap', 'i8'),
        ('proportion', 'f8'),
    ])

    profiles['pos'] = np.arange(start, stop)

    # finding depth at sites
    fragwise_overlaps = np.logical_and(
        np.greater_equal(profiles['pos'], frags['start'][:,np.newaxis]),
        np.less(profiles['pos'], frags['stop'][:,np.newaxis])
    )
    profiles['overlaps'] = np.sum(fragwise_overlaps, axis=0)

    # TODO: finding number of 5' ends at sites

    return NotImplementedError