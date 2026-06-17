"""
Genome annotation helpers (UCSC gap tracks).
"""
from .gaps import (
    GenomeGaps,
    ContigGaps,
    ucsc_hg19_gap_bed,
    b37_gap_bed,
    ucsc_hg38_gap_bed,
)

__all__ = [
    "GenomeGaps",
    "ContigGaps",
    "ucsc_hg19_gap_bed",
    "b37_gap_bed",
    "ucsc_hg38_gap_bed",
]
