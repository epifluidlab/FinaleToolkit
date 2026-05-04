from __future__ import annotations

import os
from pathlib import Path
from sys import stderr
from typing import Generator

from finaletoolkit.utils.typing import FragFile

from ._comparison import _none_geq, _none_leq
from ..io.alignment import AlignmentWrapper


def frag_generator(
    input_file: FragFile,
    contig: str | None,
    quality_threshold: int = 30,
    start: int | None = None,
    stop: int | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    intersect_policy: str = "midpoint",
    verbose: bool | int = False,
    reference_file: str | Path | None = None
) -> Generator[tuple]:
    """
    Reads from BAM, CRAM, Fragment file and returns tuples containing
    contig (chromosome), start, stop (end), mapq, and strand for each fragment.
    Optionally may filter for mapq, size, and intersection with a region.

    Parameters
    ----------
    input_file : str, pathlike, pysam TabixFile, or pysam AlignmentFile
        Fragment coordinates stored as a BAM, CRAM, or tabix-indexed
        FinaleDB fragment file. Can also be a pysam object of these files.
    contig : str or None
        Chromosome to fetch fragments over. May be None for genome-wide.
    quality_threshold : int, optional
    start : int, optional
        Left-most coordinate of interval to fetch from. See intersect_policy.
    stop : int, optional
        Right-most coordinate of interval to fetch from. See intersect_policy.
    min_length : int, optional
        Specifies lowest fragment length included in array.
    max_length : int, optional
        Specifies highest fragment length included in array.
    intersect_policy : str, optional
        Specifies what policy is used to include fragments in the
        given interval. Default is "midpoint". Policies include:
        - midpoint: the average of end coordinates of a fragment lies
        in the interval.
        - any: any part of the fragment is in the interval.
    verbose : bool, optional
    reference_file : str or Path, optional
        Reference genome file (required for CRAM).

    Returns
    -------
    frag_ends : Generator
        Generator that yields tuples:
        (contig: str, read_start: int, read_stop: int, mapq: int,
        read_on_plus: boolean)
    """
    # setting filter based on intersect policy
    if intersect_policy == 'midpoint':
        def check_intersect(r_start, r_stop, f_start, f_stop):
            return (
                (r_start is None or _none_geq(((f_start + f_stop) // 2), r_start))
                and (r_stop is None or ((f_start + f_stop) // 2) < r_stop)
            )
    elif intersect_policy == 'any':
        def check_intersect(r_start, r_stop, f_start, f_stop):
            return (
                (r_start is None or f_stop > r_start)
                and (r_stop is None or f_start < r_stop)
            )
    else:
        raise ValueError(f'{intersect_policy} is not a valid policy')

    # Raise exception if start and stop specified but not contig
    if contig is None and not (start is None and stop is None):
        if contig is None and start == 0 and stop is None:
            pass
        else:
            raise ValueError("contig should be specified if start or stop given.")

    # Use AlignmentWrapper to handle file-type-specific logic
    with AlignmentWrapper(
        input_file, 
        reference_file=reference_file, 
        quality_threshold=quality_threshold
    ) as wrapper:
        for frag in wrapper.fetch(contig, start, stop):
            try:
                if (
                    _none_geq(frag.length, min_length)
                    and _none_leq(frag.length, max_length)
                    and check_intersect(start, stop, frag.start, frag.stop)
                ):
                    yield frag.contig, frag.start, frag.stop, frag.mapq, frag.is_forward
            except TypeError as e:
                stderr.writelines([
                    "Type error encountered.\n",
                    f"Fragment length: {frag.length}\n",
                    f"fraction_low: {min_length}\n",
                    f"fraction_high: {max_length}\n",
                    "Skipping interval.\n",
                    f"Error: {e}\n"
                ])
