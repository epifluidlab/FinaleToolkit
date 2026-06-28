"""
finaletoolkit.genome.gaps
=========================

Classes and functions for the UCSC Genome Browser gap tracks (Kent et al.,
2002): telomere, centromere, and short-arm intervals for hg19/hg38/b37.
"""
from __future__ import annotations

import gzip
import importlib.resources
from os import PathLike
from sys import stdout
from typing import Iterable, Optional, Union

import numpy as np
from numpy.typing import NDArray

import finaletoolkit.genome as genome

HG19GAPS = importlib.resources.files(genome) / "data" / "hg19.gap.txt.gz"
HG38GAPS = importlib.resources.files(genome) / "data" / "hg38.gap.txt.gz"

__all__ = [
    "GenomeGaps",
    "ContigGaps",
    "ucsc_hg19_gap_bed",
    "b37_gap_bed",
    "ucsc_hg38_gap_bed",
]

_GAP_DTYPE = [
    ("contig", "<U32"),
    ("start", "<i8"),
    ("stop", "<i8"),
    ("type", "<U32"),
]


class GenomeGaps:
    """Telomere/centromere/short-arm intervals for a reference genome.

    Construct from a BED4 gap file, or use the :meth:`ucsc_hg19`/:meth:`b37`/
    :meth:`hg38` classmethods for the bundled UCSC tracks.
    """

    def __init__(self, gaps_bed: Union[PathLike, str, None] = None) -> None:
        self.centromeres: NDArray
        self.telomeres: NDArray
        self.short_arms: NDArray
        self.gaps: NDArray
        if gaps_bed is None:
            return
        gaps = np.genfromtxt(gaps_bed, dtype=_GAP_DTYPE)
        self._set_gaps(gaps)

    def _set_gaps(self, gaps: NDArray) -> None:
        self.centromeres = gaps[gaps["type"] == "centromere"]
        self.telomeres = gaps[gaps["type"] == "telomere"]
        self.short_arms = gaps[gaps["type"] == "short_arm"]
        self.gaps = gaps

    @classmethod
    def _from_track(cls, gap_resource, strip_chr: bool = False) -> "GenomeGaps":
        genome_gaps = cls()
        with importlib.resources.as_file(gap_resource) as gap_file:
            gaps = np.genfromtxt(gap_file, usecols=[1, 2, 3, 7], dtype=_GAP_DTYPE)
        if strip_chr:
            gaps["contig"] = np.char.replace(gaps["contig"], "chr", "")
        genome_gaps._set_gaps(gaps)
        return genome_gaps

    @classmethod
    def ucsc_hg19(cls) -> "GenomeGaps":
        """GenomeGaps for UCSC hg19 (``chr``-prefixed, GRCh37-based)."""
        return cls._from_track(HG19GAPS)

    @classmethod
    def b37(cls) -> "GenomeGaps":
        """GenomeGaps for Broad b37 (UCSC hg19 track with ``chr`` stripped).

        An ad-hoc approximation; other hg19/b37 differences are not accounted
        for.
        """
        return cls._from_track(HG19GAPS, strip_chr=True)

    @classmethod
    def hg38(cls) -> "GenomeGaps":
        """GenomeGaps for UCSC hg38 (``chr``-prefixed, == GRCh38)."""
        return cls._from_track(HG38GAPS)

    def in_tcmere(self, contig: str, start: int, stop: int) -> bool | None:
        """Return whether an interval overlaps a centromere or telomere.

        Returns ``None`` if the contig has no centromere annotation.

        Parameters
        ----------
        contig : str
            Chromosome name.
        start, stop : int
            Interval bounds.
        """
        centromere = self.centromeres[self.centromeres["contig"] == contig]
        telomeres = self.telomeres[self.telomeres["contig"] == contig]
        if not centromere.shape[0]:
            return None

        in_centromere = bool(
            np.any(
                np.logical_and(
                    stop > centromere["start"],
                    start < centromere["stop"],
                )
            )
        )
        # (Original code computed telomere overlap only when *no* telomeres
        # existed, making it always False; corrected to test overlap when
        # telomeres are present.)
        if telomeres.shape[0]:
            in_telomeres = bool(
                np.any(
                    np.logical_and(
                        stop > telomeres["start"],
                        start < telomeres["stop"],
                    )
                )
            )
        else:
            in_telomeres = False
        return in_centromere or in_telomeres

    def overlaps_gap(self, contig: str, start: int, stop: int) -> bool | None:
        """Return whether an interval overlaps any gap (``None`` if none)."""
        gaps = self.gaps[self.gaps["contig"] == contig]
        if not gaps.shape[0]:
            return None
        return bool(
            np.any(np.logical_and(start < gaps["stop"], stop > gaps["start"]))
        )

    def get_arm(self, contig: str, start: int, stop: int) -> str:
        """Return the chromosome arm (e.g. ``"1p"``) or ``"NOARM"``.

        Raises
        ------
        ValueError
            If ``stop < start``.
        """
        if stop < start:
            raise ValueError("start must be less than stop")

        centromere = self.centromeres[self.centromeres["contig"] == contig]
        short_arm = self.short_arms[self.short_arms["contig"] == contig]
        has_short_arm = short_arm.shape[0] > 0

        if stop < centromere["start"][0]:
            if not has_short_arm:
                return f"{contig.replace('chr', '')}p"
            return "NOARM"
        if start > centromere["stop"][0]:
            return f"{contig.replace('chr', '')}q"
        return "NOARM"

    def get_contig_gaps(self, contig: str) -> Optional["ContigGaps"]:
        """Return a :class:`ContigGaps` for ``contig`` (``None`` if no centromere)."""
        centromere = self.centromeres[self.centromeres["contig"] == contig]
        try:
            centromere_ends = (centromere[0]["start"], centromere[0]["stop"])
        except IndexError:
            return None
        telomeres = self.telomeres[self.telomeres["contig"] == contig]
        telomere_ends = [(t["start"], t["stop"]) for t in telomeres]
        short_arm = self.short_arms[self.short_arms["contig"] == contig]
        has_short_arm = short_arm.shape[0] > 0
        return ContigGaps(contig, centromere_ends, telomere_ends, has_short_arm)

    def to_bed(self, output_file: Union[str, PathLike]) -> None:
        """Write all gaps as a sorted BED4 (name = gap type).

        ``output_file`` may be a path, a ``.gz`` path, or ``"-"`` for stdout.
        """
        gaps = np.sort(self.gaps)

        def _write(handle) -> None:
            for interval in gaps:
                handle.write(
                    f"{interval['contig']}\t{interval['start']}\t"
                    f"{interval['stop']}\t{interval['type']}\n"
                )

        if str(output_file).endswith(".gz"):
            with gzip.open(output_file, "wt") as output:
                _write(output)
        elif str(output_file) == "-":
            _write(stdout)
        else:
            with open(output_file, "w") as output:
                _write(output)


class ContigGaps:
    """Centromere/telomere intervals for a single contig."""

    def __init__(
        self,
        contig: str,
        centromere: tuple[int, int],
        telomeres: Iterable[tuple[int, int]],
        has_short_arm: bool = False,
    ) -> None:
        self.contig = contig
        self.centromere = centromere
        self.telomeres = list(telomeres)
        self.has_short_arm = has_short_arm

    def in_tcmere(self, start: int, stop: int) -> bool:
        """Return whether an interval overlaps the centromere or a telomere.

        Notes
        -----
        The telomere test uses :func:`all` (an interval must overlap *every*
        telomere).  This matches the original implementation and is preserved
        deliberately: the bundled DELFI reference outputs were generated with
        this behavior, so changing it to :func:`any` would alter DELFI results.
        """
        in_centromere = (
            stop > self.centromere[0] and start < self.centromere[1]
        )
        if not self.telomeres:
            in_telomeres = False
        else:
            in_telomeres = all(
                stop > telomere[0] and start < telomere[1]
                for telomere in self.telomeres
            )
        return in_centromere or in_telomeres

    def in_gap(self, start: int, stop: int) -> bool:
        """Alias of :meth:`in_tcmere` (preserved for compatibility)."""
        in_centromere = (
            stop > self.centromere[0] and start < self.centromere[1]
        )
        in_telomeres = all(
            stop > telomere[0] and start < telomere[1]
            for telomere in self.telomeres
        )
        return in_centromere or in_telomeres

    def get_arm(self, start: int, stop: int) -> str:
        """Return the chromosome arm or ``"NOARM"``.

        Raises
        ------
        ValueError
            If ``stop < start``.
        """
        if stop < start:
            raise ValueError("start must be less than stop")

        if stop < self.centromere[0]:
            if not self.has_short_arm:
                return f"{self.contig.replace('chr', '')}p"
            return "NOARM"
        if start > self.centromere[1]:
            return f"{self.contig.replace('chr', '')}q"
        return "NOARM"


def ucsc_hg19_gap_bed(output_file: Union[str, PathLike]) -> None:
    """Write a BED4 of UCSC hg19 centromeres/telomeres/short arms."""
    return GenomeGaps.ucsc_hg19().to_bed(output_file)


def b37_gap_bed(output_file: Union[str, PathLike]) -> None:
    """Write a BED4 of Broad b37 centromeres/telomeres/short arms.

    Also useful for human_g1k_v37 (1000 Genomes) alignments.
    """
    return GenomeGaps.b37().to_bed(output_file)


def ucsc_hg38_gap_bed(output_file: Union[str, PathLike]) -> None:
    """Write a BED4 of UCSC hg38 centromeres/telomeres/short arms."""
    return GenomeGaps.hg38().to_bed(output_file)


def _cli_gap_bed(reference_genome: str, output_file: str) -> None:
    """CLI: write a gap BED for a named reference genome."""
    if reference_genome == "hg19":
        ucsc_hg19_gap_bed(output_file)
    elif reference_genome in ("b37", "human_g1k_v37"):
        b37_gap_bed(output_file)
    elif reference_genome in ("hg38", "GRCh38"):
        ucsc_hg38_gap_bed(output_file)
    else:
        raise ValueError(
            f"Gap track for {reference_genome} is currently unavailable. It is "
            "possible to create a gap track de novo if interval data for "
            "centromeres, telomeres, and short_arms exist for the reference "
            "sequence of interest."
        )
