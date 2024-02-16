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
    fraction_low: int=1,
    fraction_high: int=10000000,
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

    positions = np.arange(start, stop)

    # finding depth at sites
    fragwise_overlaps = np.logical_and(
        np.greater_equal(positions[np.newaxis], frags['start'][:,np.newaxis]),
        np.less(positions[np.newaxis], frags['stop'][:,np.newaxis])
    )
    depth = np.sum(fragwise_overlaps, axis=0)

    # finding ends
    forward_ends = np.logical_and(
        np.equal(
            positions[np.newaxis], frags['start'][:, np.newaxis]
        ), frags['strand'][:, np.newaxis]
    )
    reverse_ends = np.logical_and(
        np.equal(
            positions[np.newaxis], frags['stop'][:, np.newaxis]
        ), np.logical_not(frags['strand'][:, np.newaxis])
    )
    ends = np.sum(np.logical_or(forward_ends, reverse_ends), axis=0)
    proportions = ends/depth*100

    results = np.zeros_like(proportions, dtype=[
        ('pos', 'i8'),
        ('proportion', 'f8'),
    ])

    return proportions