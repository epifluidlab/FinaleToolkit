from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Generator, NamedTuple, Optional
import warnings

import pysam
from ..utils.logging import get_logger

logger = get_logger(__name__)

class Fragment(NamedTuple):
    """
    A zero-cost, standardized representation of a genomic fragment.
    This NamedTuple provides a consistent interface regardless of whether
    the data originated from a BAM/CRAM record or a BED/Tabix line.
    """
    contig: str
    start: int
    stop: int
    mapq: int
    is_forward: bool

    @property
    def length(self) -> int:
        """Returns the length of the fragment."""
        return self.stop - self.start

class AlignmentWrapper:
    """
    A unified interface for reading genomic alignment and fragment data.
    Wraps BAM, CRAM, and tabix-indexed fragment files (frag.gz).
    
    This class provides a consistent way to access contig information and
    iterate over fragments across different file formats.
    """

    def __init__(
        self, 
        path: str | Path | pysam.AlignmentFile | pysam.TabixFile,
        reference_file: Optional[str | Path] = None,
        threads: int = 1,
        quality_threshold: int = 30,
        read1_only: bool = True
    ):
        """
        Initialize the AlignmentWrapper.

        Parameters
        ----------
        path : str | Path | pysam.AlignmentFile | pysam.TabixFile
            Path to the BAM, CRAM, or Frag.gz file, or an already-open
            pysam alignment or tabix handle.
        reference_file : str | Path, optional
            Path to the reference genome (required for CRAM).
        threads : int, optional
            Number of threads to use for decompression (BAM/CRAM).
        quality_threshold : int, optional
            Minimum mapping quality for fragments. Default is 30.
        read1_only : bool, optional
            If True, only processes read1 from BAM/CRAM files to avoid 
            double-counting fragments. Default is True.
        """
        self.path = str(path) if isinstance(path, (str, Path)) else None
        self.reference_file = str(reference_file) if reference_file else None
        self.threads = threads
        self.quality_threshold = quality_threshold
        self.read1_only = read1_only
        self._handle = None
        self._is_sam = False
        self._is_tabix = False
        self._chroms = None
        self._bed_format = False # Used for tabix files
        self._owns_handle = False

        if isinstance(path, pysam.AlignmentFile):
            self._handle = path
            self._is_sam = True
            self._chroms = dict(zip(self._handle.references, self._handle.lengths))
            return
        if isinstance(path, pysam.TabixFile):
            self._handle = path
            self._is_tabix = True
            self._chroms = {c: None for c in self._handle.contigs}
            self._detect_tabix_format(self._handle)
            return

        if self.path is None or not os.path.exists(self.path):
            raise FileNotFoundError(f"Alignment file not found: {path}")

        self._open_file()

    def _detect_tabix_format(self, tabix_handle: pysam.TabixFile):
        try:
            first_line = next(tabix_handle.fetch(parser=pysam.asTuple()))
            if len(first_line) > 5:
                self._bed_format = True
                warnings.warn(
                    "input_file does not follow Fragmentation file format "
                    "accepted by FinaleToolkit. Attempting to read as a BED6 "
                    "file.",
                    UserWarning
                )
        except StopIteration:
            pass

    def _open_file(self):
        """Detects file type and opens the appropriate handle."""
        lower_path = self.path.lower()
        
        if lower_path.endswith((".bam", ".cram", ".sam")):
            self._is_sam = True
            
            # Check for index if not a SAM file
            if lower_path.endswith(".bam"):
                if not (os.path.exists(self.path + ".bai") or os.path.exists(self.path[:-4] + ".bai")):
                    raise FileNotFoundError(f"BAM file {self.path} missing index (.bai)")
            elif lower_path.endswith(".cram"):
                if not (os.path.exists(self.path + ".crai") or os.path.exists(self.path[:-5] + ".crai")):
                    raise FileNotFoundError(f"CRAM file {self.path} missing index (.crai)")

            # For CRAM, reference_filename is often necessary
            self._handle = pysam.AlignmentFile(
                self.path, "r", 
                reference_filename=self.reference_file,
                threads=self.threads
            )
            self._owns_handle = True
            self._chroms = dict(zip(self._handle.references, self._handle.lengths))
            
        elif lower_path.endswith((".gz", ".bgz")):
            # Check for tabix index
            if os.path.exists(self.path + ".tbi"):
                self._is_tabix = True
                self._handle = pysam.TabixFile(self.path)
                self._owns_handle = True
                self._chroms = {c: None for c in self._handle.contigs}
                self._detect_tabix_format(self._handle)
            else:
                raise FileNotFoundError(f"Compressed file {self.path} missing tabix index (.tbi)")
        else:
            raise ValueError(f"Unsupported file format: {self.path}")

    @property
    def chroms(self) -> Dict[str, Optional[int]]:
        """Returns a dictionary of contig names and their lengths (if available)."""
        return self._chroms

    @property
    def is_sam(self) -> bool:
        """Returns True if the source is a SAM/BAM/CRAM file."""
        return self._is_sam

    def fetch(
        self, 
        contig: Optional[str] = None, 
        start: Optional[int] = None, 
        stop: Optional[int] = None
    ) -> Generator[Fragment, None, None]:
        """
        Unified iterator over standardized Fragment objects in a region.
        """
        if self._is_sam:
            for read in self._handle.fetch(contig, start, stop):
                # 1. Basic Quality Filtering (Equivalent to low_quality_read_pairs)
                if (read.is_unmapped or 
                    read.is_secondary or 
                    not read.is_paired or 
                    read.mate_is_unmapped or 
                    read.is_duplicate or 
                    read.mapping_quality < self.quality_threshold or 
                    read.is_qcfail or 
                    read.is_supplementary or 
                    not read.is_proper_pair):
                    continue
                
                # 2. Avoid double-counting (Read1 only)
                if self.read1_only and read.is_read2:
                    continue

                # 3. Calculate fragment span
                if read.template_length > 0:
                    f_start = read.reference_start
                    f_stop = read.reference_start + read.template_length
                elif read.template_length < 0:
                    f_start = read.reference_end + read.template_length
                    f_stop = read.reference_end
                else:
                    continue # Skip if template length is 0

                yield Fragment(
                    contig=read.reference_name,
                    start=f_start,
                    stop=f_stop,
                    mapq=read.mapping_quality,
                    is_forward=read.is_forward
                )

        elif self._is_tabix:
            for line in self._handle.fetch(contig, start, stop, parser=pysam.asTuple()):
                try:
                    f_start = int(line[1])
                    f_stop = int(line[2])
                    
                    if self._bed_format:
                        mapq = int(line[4])
                        is_forward = '+' in line[5]
                    else:
                        mapq = int(line[3])
                        is_forward = '+' in line[4]

                    if mapq < self.quality_threshold:
                        continue

                    yield Fragment(
                        contig=line[0],
                        start=f_start,
                        stop=f_stop,
                        mapq=mapq,
                        is_forward=is_forward
                    )
                except (ValueError, IndexError):
                    continue

    def close(self):
        """Closes the underlying file handle."""
        if self._handle and self._owns_handle:
            self._handle.close()
            self._handle = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()
