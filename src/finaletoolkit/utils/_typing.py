"""
# _typing
Some useful type aliases
"""
from __future__ import annotations
from os import PathLike
from typing import Union

from pysam import AlignmentFile, TabixFile

# files accepted by frag_generator
FragFile = Union[str, PathLike, AlignmentFile, TabixFile]
Intervals = Union[str, PathLike]
