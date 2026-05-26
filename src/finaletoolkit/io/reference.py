from __future__ import annotations

import os
import threading
from pathlib import Path
from typing import Dict

import py2bit
import pysam

from finaletoolkit.utils.logging import get_logger

logger = get_logger(__name__)

class ReferenceWrapper:
    """
    Wraps a 2-bit or FASTA-formatted reference for queries. This encapsulates most of the
    py2bit/pysam logic behind a common interface.
    
    This class is optionally thread-safe via a lock. The class can be used within a context manager to handle scope/lifetime.
    """

    def __init__(self, reference_path: str | Path, use_lock: bool = True):
        """
        Initialize the ReferenceWrapper with a path to a 2-bit or FASTA file.

        Parameters
        ----------
        reference_path : str | Path
            Path to the reference file.
        use_lock : bool, optional
            If True, uses a threading lock for sequence queries. This is necessary
            if the same instance is shared across threads. If False, no lock is used,
            which is faster but requires the user to manage thread-safety (e.g., by
            giving each thread its own instance or using thread-local storage).
            Defaults to True.
        """
        self.reference_path = str(reference_path)
        self.use_lock = use_lock
        self._lock = threading.Lock() if use_lock else None
        self._handle = None
        self._is_2bit = False
        self._is_fasta = False
        self._chroms = None

        if not os.path.exists(self.reference_path):
            error_msg = f"Reference file not found: {self.reference_path}"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)

        if self.reference_path.endswith((".2bit", ".tb2")):
            self._is_2bit = True
            self._open_2bit()
        elif self.reference_path.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")):
            self._is_fasta = True
            self._open_fasta()
        else:
            # Try to infer or default to fasta if unknown but exists
            logger.warning(f"Unknown reference file extension for {self.reference_path}. Attempting to open as FASTA.")
            self._is_fasta = True
            self._open_fasta()

    def _open_2bit(self):
        try:
            self._handle = py2bit.open(self.reference_path)
            self._chroms = self._handle.chroms()
        except Exception as e:
            logger.error(f"Failed to open 2-bit file {self.reference_path}: {e}")
            raise

    def _open_fasta(self):
        fai_path = self.reference_path + ".fai"
        if not os.path.exists(fai_path):
            logger.info(
                f"FASTA index not found for {self.reference_path}. "
                "Creating index with pysam.faidx()..."
            )
            pysam.faidx(self.reference_path)
        try:
            self._handle = pysam.FastaFile(self.reference_path)
            self._chroms = dict(zip(self._handle.references, self._handle.lengths))
        except Exception as e:
            logger.error(f"Failed to open FASTA file {self.reference_path}: {e}")
            raise

    @property
    def chroms(self) -> Dict[str, int]:
        """
        Returns a dictionary mapping chromosome names to their lengths.
        """
        return self._chroms

    def sequence(self, contig: str,
                 start: int | None = None,
                 stop: int | None = None,
                 fail_on_excess_range: bool = True
                 ) -> str:
        """
        Returns the sequence for the queried region / contig.
        
        Parameters
        ----------
        contig : str
            Chromosome or contig name.
        start : int, optional
            0-based start position.
        stop : int, optional
            0-based end position (exclusive).
        fail_on_excess_range : bool, optional
            If True, raises ValueError if the requested range is outside chromosome bounds.
            If False, truncates the range to chromosome bounds.
        """
        if contig not in self._chroms:
            raise ValueError(f"Contig {contig} not found in reference.")

        chrom_len = self._chroms[contig]
        
        if start is None:
            start = 0
        if stop is None:
            stop = chrom_len

        if start < 0 or stop > chrom_len or start > stop:
            if fail_on_excess_range:
                raise ValueError(f"Requested range {contig}:{start}-{stop} is out of bounds (0-{chrom_len}).")
            else:
                start = max(0, start)
                stop = min(chrom_len, stop)
                if start > stop:
                    return ""

        if self.use_lock:
            with self._lock:
                if self._is_2bit:
                    # py2bit.sequence(chrom, start, end) is 0-based, half-open
                    seq = self._handle.sequence(contig, start, stop)
                else:
                    # pysam.FastaFile.fetch(contig, start, stop) is 0-based, half-open
                    seq = self._handle.fetch(contig, start, stop)
        else:
            if self._is_2bit:
                seq = self._handle.sequence(contig, start, stop)
            else:
                seq = self._handle.fetch(contig, start, stop)
        
        return seq.upper() if seq else ""

    def close(self):
        """
        Safely closes the underlying file handle.
        """
        if self.use_lock:
            with self._lock:
                if self._handle is not None:
                    self._handle.close()
                    self._handle = None
        else:
            if self._handle is not None:
                self._handle.close()
                self._handle = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __getitem__(self, contig: str) -> _ContigSlicer:
        """
        Returns a slicer object for the specified contig, allowing for
        syntax like ref['chr1'][100:200].
        """
        if contig not in self._chroms:
            raise ValueError(f"Contig {contig} not found in reference.")
        return _ContigSlicer(self, contig)

class _ContigSlicer:
    """
    Helper class to enable slicing syntax for a specific contig.
    """
    def __init__(self, wrapper: ReferenceWrapper, contig: str):
        self.wrapper = wrapper
        self.contig = contig

    def __getitem__(self, key: slice | int) -> str:
        if isinstance(key, slice):
            return self.wrapper.sequence(self.contig, key.start, key.stop)
        elif isinstance(key, int):
            return self.wrapper.sequence(self.contig, key, key + 1)
        else:
            raise TypeError("Slicer indices must be integers or slices.")

    def __len__(self) -> int:
        return self.wrapper.chroms[self.contig]
