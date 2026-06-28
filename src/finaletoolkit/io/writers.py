"""
Small, shared output helpers used by the feature writers.

These centralize the "open a path that might be ``-`` (stdout) or ``.gz``"
boilerplate that was duplicated across the original feature modules, while
leaving the exact line formatting to each feature so output stays byte-for-byte
compatible with the original toolkit.
"""
from __future__ import annotations

import gzip
import sys
from contextlib import contextmanager
from typing import IO, Iterator

__all__ = ["smart_open_text", "is_stdout"]


def is_stdout(path: str) -> bool:
    """Return ``True`` if ``path`` denotes standard output (``"-"``)."""
    return path == "-"


@contextmanager
def smart_open_text(path: str, mode: str = "wt") -> Iterator[IO[str]]:
    """Open a text destination, transparently handling stdout and gzip.

    Parameters
    ----------
    path : str
        Output path.  ``"-"`` writes to :data:`sys.stdout`; a ``.gz`` suffix
        selects gzip compression; anything else is a plain text file.
    mode : str, optional
        Text mode, by default ``"wt"``.

    Yields
    ------
    file object
        A writable text stream.  Stdout is *not* closed on exit.
    """
    if is_stdout(path):
        yield sys.stdout
        return

    if str(path).endswith(".gz"):
        handle: IO[str] = gzip.open(path, mode)
    else:
        handle = open(path, mode)
    try:
        yield handle
    finally:
        handle.close()
