"""
Unified access to 2bit and FASTA reference genomes.

:class:`ReferenceWrapper` hides the py2bit/pysam differences behind a single
:meth:`ReferenceWrapper.sequence` method and is optionally thread-safe.
"""
from __future__ import annotations

import os
import threading
from pathlib import Path
from typing import Dict

import py2bit
import pysam

from ..exceptions import ContigNotFoundError, OutOfBoundsError
from ..utils.logging import get_logger

logger = get_logger(__name__)

__all__ = ["ReferenceWrapper"]

_TWOBIT_SUFFIXES = (".2bit", ".tb2")
_FASTA_SUFFIXES = (
    ".fa",
    ".fasta",
    ".fna",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
)


class ReferenceWrapper:
    """Wrap a 2bit- or FASTA-formatted reference genome for sequence queries.

    The wrapper encapsulates the py2bit/pysam handle logic behind a common
    interface and may be used as a context manager.

    Parameters
    ----------
    reference_path : str or Path
        Path to a ``.2bit``/``.tb2`` or FASTA (``.fa``/``.fasta``/``.fna``,
        optionally ``.gz``) reference file.
    use_lock : bool, optional
        If ``True`` (default), guard sequence queries with a
        :class:`threading.Lock` so a single instance may be shared across
        threads.  Set ``False`` for per-thread/per-process instances to avoid
        lock overhead.

    Raises
    ------
    FileNotFoundError
        If ``reference_path`` does not exist.
    """

    def __init__(self, reference_path: str | Path, use_lock: bool = True) -> None:
        self.reference_path = str(reference_path)
        self.use_lock = use_lock
        self._lock = threading.Lock() if use_lock else None
        self._handle = None
        self._is_2bit = False
        self._is_fasta = False
        self._chroms: Dict[str, int] | None = None

        if not os.path.exists(self.reference_path):
            error_msg = f"Reference file not found: {self.reference_path}"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)

        if self.reference_path.endswith(_TWOBIT_SUFFIXES):
            self._is_2bit = True
            self._open_2bit()
        elif self.reference_path.endswith(_FASTA_SUFFIXES):
            self._is_fasta = True
            self._open_fasta()
        else:
            logger.warning(
                "Unknown reference file extension for %s. "
                "Attempting to open as FASTA.",
                self.reference_path,
            )
            self._is_fasta = True
            self._open_fasta()

    # -- handle management --------------------------------------------------

    def _open_2bit(self) -> None:
        try:
            self._handle = py2bit.open(self.reference_path)
            self._chroms = self._handle.chroms()
        except Exception as e:
            logger.error("Failed to open 2-bit file %s: %s", self.reference_path, e)
            raise

    def _open_fasta(self) -> None:
        fai_path = self.reference_path + ".fai"
        if not os.path.exists(fai_path):
            logger.info(
                "FASTA index not found for %s. Creating index with "
                "pysam.faidx()...",
                self.reference_path,
            )
            pysam.faidx(self.reference_path)
        try:
            self._handle = pysam.FastaFile(self.reference_path)
            self._chroms = dict(zip(self._handle.references, self._handle.lengths))
        except Exception as e:
            logger.error("Failed to open FASTA file %s: %s", self.reference_path, e)
            raise

    # -- public interface ---------------------------------------------------

    @property
    def chroms(self) -> Dict[str, int]:
        """Mapping of chromosome name to length."""
        return self._chroms

    def sequence(
        self,
        contig: str,
        start: int | None = None,
        stop: int | None = None,
        fail_on_excess_range: bool = True,
    ) -> str:
        """Return the upper-cased reference sequence for a region.

        Parameters
        ----------
        contig : str
            Chromosome or contig name.
        start : int, optional
            0-based start position (defaults to 0).
        stop : int, optional
            0-based, exclusive end position (defaults to the contig length).
        fail_on_excess_range : bool, optional
            If ``True`` (default), raise when the range exceeds the contig
            bounds; otherwise truncate the range to the contig.

        Returns
        -------
        str
            The upper-cased sequence (empty string if the truncated range is
            empty).

        Raises
        ------
        ContigNotFoundError
            If ``contig`` is absent from the reference.
        OutOfBoundsError
            If the range is out of bounds and ``fail_on_excess_range`` is True.
        """
        if contig not in self._chroms:
            raise ContigNotFoundError(f"Contig {contig} not found in reference.")

        chrom_len = self._chroms[contig]

        if start is None:
            start = 0
        if stop is None:
            stop = chrom_len

        if start < 0 or stop > chrom_len or start > stop:
            if fail_on_excess_range:
                raise OutOfBoundsError(
                    f"Requested range {contig}:{start}-{stop} is out of bounds "
                    f"(0-{chrom_len})."
                )
            start = max(0, start)
            stop = min(chrom_len, stop)
            if start > stop:
                return ""

        seq = self._fetch(contig, start, stop)
        return seq.upper() if seq else ""

    def _fetch(self, contig: str, start: int, stop: int) -> str:
        """Backend-agnostic fetch, guarded by the lock when enabled."""
        if self._lock is not None:
            with self._lock:
                return self._fetch_unlocked(contig, start, stop)
        return self._fetch_unlocked(contig, start, stop)

    def _fetch_unlocked(self, contig: str, start: int, stop: int) -> str:
        # Both backends use 0-based, half-open coordinates.
        if self._is_2bit:
            return self._handle.sequence(contig, start, stop)
        return self._handle.fetch(contig, start, stop)

    def close(self) -> None:
        """Safely close the underlying file handle."""
        if self._lock is not None:
            with self._lock:
                self._close_handle()
        else:
            self._close_handle()

    def _close_handle(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    # -- context manager / slicing -----------------------------------------

    def __enter__(self) -> "ReferenceWrapper":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def __getitem__(self, contig: str) -> "_ContigSlicer":
        """Enable ``ref['chr1'][100:200]`` slicing syntax."""
        if contig not in self._chroms:
            raise ContigNotFoundError(f"Contig {contig} not found in reference.")
        return _ContigSlicer(self, contig)


class _ContigSlicer:
    """Helper enabling ``ref[contig][start:stop]`` slicing of a reference."""

    def __init__(self, wrapper: ReferenceWrapper, contig: str) -> None:
        self.wrapper = wrapper
        self.contig = contig

    def __getitem__(self, key: slice | int) -> str:
        if isinstance(key, slice):
            return self.wrapper.sequence(self.contig, key.start, key.stop)
        if isinstance(key, int):
            return self.wrapper.sequence(self.contig, key, key + 1)
        raise TypeError("Slicer indices must be integers or slices.")

    def __len__(self) -> int:
        return self.wrapper.chroms[self.contig]
