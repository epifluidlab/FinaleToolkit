"""
Informative exception hierarchy for FinaleToolkit.

Every exception below subclasses a built-in exception type (``ValueError``,
``FileNotFoundError``, ``IndexError`` ...) in addition to
:class:`FinaleToolkitError`.  This keeps the library a strict drop-in for the
original package: existing ``except ValueError``/``except FileNotFoundError``
handlers continue to catch these errors, while new code can catch the more
specific :class:`FinaleToolkitError` subclasses or test for them by type.

Examples
--------
>>> from finaletoolkit.exceptions import UnsupportedFormatError
>>> try:
...     raise UnsupportedFormatError("not a bam")
... except ValueError:  # still works -- UnsupportedFormatError is a ValueError
...     print("caught")
caught
"""
from __future__ import annotations


class FinaleToolkitError(Exception):
    """Base class for all FinaleToolkit-specific errors."""


class InvalidInputError(FinaleToolkitError, ValueError):
    """Raised when user-supplied input is malformed or inconsistent.

    Subclasses :class:`ValueError` for backwards compatibility.
    """


class UnsupportedFormatError(InvalidInputError):
    """Raised when an input file is in a format the toolkit cannot read."""


class MissingReferenceError(InvalidInputError):
    """Raised when a reference genome is required but not provided.

    For example, CRAM input requires a FASTA reference.
    """


class MissingIndexError(FinaleToolkitError, FileNotFoundError):
    """Raised when a required index (``.bai``/``.crai``/``.tbi``/``.fai``) is
    missing.  Subclasses :class:`FileNotFoundError`.
    """


class ContigNotFoundError(InvalidInputError):
    """Raised when a requested contig is absent from a reference or alignment."""


class ContigMismatchError(InvalidInputError):
    """Raised when contigs (names or sizes) of two files are incompatible."""


class OutOfBoundsError(InvalidInputError, IndexError):
    """Raised when a queried interval falls outside chromosome bounds.

    Subclasses both :class:`ValueError` (via :class:`InvalidInputError`) and
    :class:`IndexError` so that either handler catches it.
    """


__all__ = [
    "FinaleToolkitError",
    "InvalidInputError",
    "UnsupportedFormatError",
    "MissingReferenceError",
    "MissingIndexError",
    "ContigNotFoundError",
    "ContigMismatchError",
    "OutOfBoundsError",
]
