"""
I/O layer: unified reference and alignment access plus output helpers.
"""
from .reference import ReferenceWrapper
from .alignment import AlignmentWrapper, Fragment
from .writers import smart_open_text, is_stdout

__all__ = [
    "ReferenceWrapper",
    "AlignmentWrapper",
    "Fragment",
    "smart_open_text",
    "is_stdout",
]
