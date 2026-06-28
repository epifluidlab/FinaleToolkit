"""
Unified reading of BAM/CRAM/SAM and tabix-indexed fragment files.

:class:`AlignmentWrapper` yields a stream of standardized :class:`Fragment`
records regardless of source format, applying the same quality filtering used
throughout the toolkit.
"""
from __future__ import annotations

import os
import warnings
from pathlib import Path
from typing import Dict, Generator, NamedTuple, Optional

import pysam

from ..exceptions import MissingIndexError, UnsupportedFormatError
from ..utils.logging import get_logger

logger = get_logger(__name__)

__all__ = ["Fragment", "AlignmentWrapper"]


class Fragment(NamedTuple):
    """A zero-cost, standardized representation of a genomic fragment.

    The same record type is produced whether the data came from a BAM/CRAM
    record or a tabix/BED line, so downstream code never branches on format.

    Attributes
    ----------
    contig : str
        Chromosome/contig name.
    start : int
        0-based, inclusive fragment start.
    stop : int
        0-based, exclusive fragment stop.
    mapq : int
        Mapping quality.
    is_forward : bool
        ``True`` if the fragment is on the ``+`` strand.
    """

    contig: str
    start: int
    stop: int
    mapq: int
    is_forward: bool

    @property
    def length(self) -> int:
        """The fragment length in base pairs (``stop - start``)."""
        return self.stop - self.start


# samtools-style flag filtering applied to BAM/CRAM reads (equivalent to
# ``-F 3852 -f 3``): drop unmapped/secondary/duplicate/qcfail/supplementary
# reads and keep only properly-paired reads.
def _read_is_low_quality(read: "pysam.AlignedSegment", quality_threshold: int) -> bool:
    return (
        read.is_unmapped
        or read.is_secondary
        or not read.is_paired
        or read.mate_is_unmapped
        or read.is_duplicate
        or read.mapping_quality < quality_threshold
        or read.is_qcfail
        or read.is_supplementary
        or not read.is_proper_pair
    )


class AlignmentWrapper:
    """A unified interface for reading alignment and fragment data.

    Wraps BAM, CRAM, SAM, and tabix-indexed fragment files (``frag.gz`` and
    BED6 ``bed.gz``), exposing a consistent :meth:`fetch` generator and contig
    table.

    Parameters
    ----------
    path : str, Path, pysam.AlignmentFile, or pysam.TabixFile
        Path to the input file, or an already-open pysam handle.
    reference_file : str or Path, optional
        Reference genome path (required for CRAM input).
    threads : int, optional
        Decompression threads for BAM/CRAM (default 1).
    quality_threshold : int, optional
        Minimum mapping quality for emitted fragments (default 30).
    read1_only : bool, optional
        If ``True`` (default), only read1 is used from BAM/CRAM to avoid
        double-counting fragments.  Has no effect on tabix files.

    Raises
    ------
    FileNotFoundError
        If the file or a required index is missing.
    UnsupportedFormatError
        If the file extension is not recognized.
    """

    def __init__(
        self,
        path: str | Path | pysam.AlignmentFile | pysam.TabixFile,
        reference_file: Optional[str | Path] = None,
        threads: int = 1,
        quality_threshold: int = 30,
        read1_only: bool = True,
    ) -> None:
        self.path = str(path) if isinstance(path, (str, Path)) else None
        self.reference_file = str(reference_file) if reference_file else None
        self.threads = threads
        self.quality_threshold = quality_threshold
        self.read1_only = read1_only
        self._handle = None
        self._is_sam = False
        self._is_tabix = False
        self._chroms: Dict[str, Optional[int]] | None = None
        self._bed_format = False  # tabix BED6 vs FinaleDB frag layout
        self._owns_handle = False

        # Accept already-open pysam handles without taking ownership.
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

    # -- format detection ---------------------------------------------------

    def _detect_tabix_format(self, tabix_handle: pysam.TabixFile) -> None:
        """Distinguish FinaleDB fragment layout from BED6+ by column count."""
        try:
            first_line = next(tabix_handle.fetch(parser=pysam.asTuple()))
            if len(first_line) > 5:
                self._bed_format = True
                warnings.warn(
                    "input_file does not follow Fragmentation file format "
                    "accepted by FinaleToolkit. Attempting to read as a BED6 "
                    "file.",
                    UserWarning,
                )
        except StopIteration:
            pass

    def _open_file(self) -> None:
        """Detect the file type by extension and open the appropriate handle."""
        lower_path = self.path.lower()

        if lower_path.endswith((".bam", ".cram", ".sam")):
            self._is_sam = True

            if lower_path.endswith(".bam"):
                if not (
                    os.path.exists(self.path + ".bai")
                    or os.path.exists(self.path[:-4] + ".bai")
                ):
                    raise MissingIndexError(
                        f"BAM file {self.path} missing index (.bai)"
                    )
            elif lower_path.endswith(".cram"):
                if not (
                    os.path.exists(self.path + ".crai")
                    or os.path.exists(self.path[:-5] + ".crai")
                ):
                    raise MissingIndexError(
                        f"CRAM file {self.path} missing index (.crai)"
                    )

            self._handle = pysam.AlignmentFile(
                self.path,
                "r",
                reference_filename=self.reference_file,
                threads=self.threads,
            )
            self._owns_handle = True
            self._chroms = dict(zip(self._handle.references, self._handle.lengths))

        elif lower_path.endswith((".gz", ".bgz")):
            if os.path.exists(self.path + ".tbi"):
                self._is_tabix = True
                self._handle = pysam.TabixFile(self.path)
                self._owns_handle = True
                self._chroms = {c: None for c in self._handle.contigs}
                self._detect_tabix_format(self._handle)
            else:
                raise MissingIndexError(
                    f"Compressed file {self.path} missing tabix index (.tbi)"
                )
        else:
            raise UnsupportedFormatError(f"Unsupported file format: {self.path}")

    # -- public interface ---------------------------------------------------

    @property
    def chroms(self) -> Dict[str, Optional[int]]:
        """Mapping of contig name to length (``None`` for tabix files)."""
        return self._chroms

    @property
    def is_sam(self) -> bool:
        """``True`` if the source is a SAM/BAM/CRAM file."""
        return self._is_sam

    def fetch(
        self,
        contig: Optional[str] = None,
        start: Optional[int] = None,
        stop: Optional[int] = None,
    ) -> Generator[Fragment, None, None]:
        """Yield :class:`Fragment` records overlapping a region.

        Parameters
        ----------
        contig : str, optional
            Contig to restrict to (``None`` for the whole file).
        start, stop : int, optional
            Region bounds passed to the underlying index query.

        Yields
        ------
        Fragment
            Standardized fragment records passing the quality filter.
        """
        if self._is_sam:
            yield from self._fetch_sam(contig, start, stop)
        elif self._is_tabix:
            yield from self._fetch_tabix(contig, start, stop)

    def _fetch_sam(self, contig, start, stop) -> Generator[Fragment, None, None]:
        read1_only = self.read1_only
        quality_threshold = self.quality_threshold
        for read in self._handle.fetch(contig, start, stop):
            if _read_is_low_quality(read, quality_threshold):
                continue
            if read1_only and read.is_read2:
                continue

            # Reconstruct the fragment span from the template length.
            tlen = read.template_length
            if tlen > 0:
                f_start = read.reference_start
                f_stop = read.reference_start + tlen
            elif tlen < 0:
                f_start = read.reference_end + tlen
                f_stop = read.reference_end
            else:
                continue

            yield Fragment(
                contig=read.reference_name,
                start=f_start,
                stop=f_stop,
                mapq=read.mapping_quality,
                is_forward=read.is_forward,
            )

    def _fetch_tabix(self, contig, start, stop) -> Generator[Fragment, None, None]:
        quality_threshold = self.quality_threshold
        bed_format = self._bed_format
        for line in self._handle.fetch(
            contig,
            start,
            stop,
            parser=pysam.asTuple(),
            multiple_iterators=True,
        ):
            try:
                f_start = int(line[1])
                f_stop = int(line[2])

                if bed_format:
                    mapq = int(line[4])
                    is_forward = "+" in line[5]
                else:
                    mapq = int(line[3])
                    is_forward = "+" in line[4]

                if mapq < quality_threshold:
                    continue

                yield Fragment(
                    contig=line[0],
                    start=f_start,
                    stop=f_stop,
                    mapq=mapq,
                    is_forward=is_forward,
                )
            except (ValueError, IndexError):
                continue

    def close(self) -> None:
        """Close the underlying handle if this wrapper owns it."""
        if self._handle and self._owns_handle:
            self._handle.close()
            self._handle = None

    def __enter__(self) -> "AlignmentWrapper":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass
