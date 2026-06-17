"""
Shared type aliases used throughout FinaleToolkit.
"""
from __future__ import annotations

from os import PathLike
from typing import Union

from pysam import AlignmentFile, TabixFile

# Files accepted by ``frag_generator`` / ``frag_array``: a path to a BAM, CRAM,
# or tabix-indexed fragment file, or an already-open pysam handle.
FragFile = Union[str, PathLike, AlignmentFile, TabixFile]

# A ``.chrom.sizes`` file (tab-delimited contig name / length).
ChromSizes = Union[str, PathLike]

# A BED file of intervals.
Intervals = Union[str, PathLike]

__all__ = ["FragFile", "ChromSizes", "Intervals"]
