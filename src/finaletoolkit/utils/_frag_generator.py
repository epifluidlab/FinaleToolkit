"""
Streaming fragment generator shared by every fragmentation feature.
"""
from __future__ import annotations

from pathlib import Path
from sys import stderr
from typing import Callable, Generator, Tuple

from ..exceptions import InvalidInputError
from ..io.alignment import AlignmentWrapper
from ._comparison import _none_geq, _none_leq
from .typing import FragFile

__all__ = ["frag_generator"]

# Yielded record: (contig, start, stop, mapq, is_forward)
FragTuple = Tuple[str, int, int, int, bool]


def _make_intersect_checker(intersect_policy: str) -> Callable[[int | None, int | None, int, int], bool]:
    """Return the region-membership predicate for an intersect policy.

    Parameters
    ----------
    intersect_policy : {"midpoint", "any"}
        ``"midpoint"`` includes a fragment when its midpoint lies in
        ``[start, stop)``; ``"any"`` includes it when any base overlaps.

    Raises
    ------
    InvalidInputError
        If the policy name is not recognized.
    """
    if intersect_policy == "midpoint":

        def check_intersect(r_start, r_stop, f_start, f_stop) -> bool:
            midpoint = (f_start + f_stop) // 2
            return (
                (r_start is None or _none_geq(midpoint, r_start))
                and (r_stop is None or midpoint < r_stop)
            )

    elif intersect_policy == "any":

        def check_intersect(r_start, r_stop, f_start, f_stop) -> bool:
            return (
                (r_start is None or f_stop > r_start)
                and (r_stop is None or f_start < r_stop)
            )

    else:
        raise InvalidInputError(f"{intersect_policy} is not a valid policy")

    return check_intersect


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
    reference_file: str | Path | None = None,
) -> Generator[FragTuple, None, None]:
    """Stream fragments from a BAM/CRAM/fragment file with optional filtering.

    Parameters
    ----------
    input_file : str, Path, pysam.TabixFile, or pysam.AlignmentFile
        Fragment coordinates as a BAM, CRAM, or tabix-indexed FinaleDB
        fragment file (or an open pysam handle).
    contig : str or None
        Chromosome to fetch over (``None`` for genome-wide).
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    start, stop : int, optional
        Region bounds; see ``intersect_policy``.
    min_length, max_length : int, optional
        Inclusive fragment-length bounds (``None`` disables a bound).
    intersect_policy : {"midpoint", "any"}, optional
        Region-membership policy (default ``"midpoint"``).
    verbose : bool or int, optional
        Unused placeholder kept for signature compatibility.
    reference_file : str or Path, optional
        Reference genome (required for CRAM).

    Yields
    ------
    tuple
        ``(contig, start, stop, mapq, is_forward)`` for each passing fragment.

    Raises
    ------
    InvalidInputError
        If ``start``/``stop`` are given without ``contig`` (except the
        whole-genome special case ``start == 0, stop is None``).
    """
    check_intersect = _make_intersect_checker(intersect_policy)

    # Require a contig when bounds are given (whole-genome start==0/stop==None ok).
    if contig is None and not (start is None and stop is None):
        if not (start == 0 and stop is None):
            raise InvalidInputError(
                "contig should be specified if start or stop given."
            )

    with AlignmentWrapper(
        input_file,
        reference_file=reference_file,
        quality_threshold=quality_threshold,
    ) as wrapper:
        for frag in wrapper.fetch(contig, start, stop):
            try:
                if (
                    _none_geq(frag.length, min_length)
                    and _none_leq(frag.length, max_length)
                    and check_intersect(start, stop, frag.start, frag.stop)
                ):
                    yield (
                        frag.contig,
                        frag.start,
                        frag.stop,
                        frag.mapq,
                        frag.is_forward,
                    )
            except TypeError as e:
                stderr.writelines(
                    [
                        "Type error encountered.\n",
                        f"Fragment length: {frag.length}\n",
                        f"fraction_low: {min_length}\n",
                        f"fraction_high: {max_length}\n",
                        "Skipping interval.\n",
                        f"Error: {e}\n",
                    ]
                )
