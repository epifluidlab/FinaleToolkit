"""
FinaleToolkit
=============

A package and standalone program to extract fragmentation features of
cell-free DNA from paired-end sequencing data.

This rewrite preserves the original public API exactly while adding a single,
flat import namespace so that every feature can be reached without knowing the
internal submodule layout::

    import finaletoolkit as ftk

    cov = ftk.coverage("sample.bam", "intervals.bed", "-")
    motifs = ftk.end_motifs("sample.bam", "hg38.2bit")
    gaps = ftk.GenomeGaps.hg38()

The legacy subpackages remain importable and unchanged::

    from finaletoolkit.frag import wps, multi_wps
    from finaletoolkit.utils import frag_generator
    from finaletoolkit.genome import GenomeGaps

Heavy dependencies (pysam, numba, pandas ...) are imported lazily on first
access so that ``import finaletoolkit`` stays cheap.
"""
from __future__ import annotations

from .version import __version__
from . import exceptions
from .exceptions import (
    ContigMismatchError,
    ContigNotFoundError,
    FinaleToolkitError,
    InvalidInputError,
    MissingIndexError,
    MissingReferenceError,
    OutOfBoundsError,
    UnsupportedFormatError,
)

# Subpackages that may be accessed as attributes (lazy-loaded, as in the
# original package).
_SUBMODULES = ("cli", "frag", "genome", "io", "utils")

# Flat namespace: maps each public symbol to the submodule that defines it.
# ``finaletoolkit.<name>`` lazily imports ``finaletoolkit.<submodule>`` and
# returns ``getattr(submodule, <attr>)``.
_EXPORTS: dict[str, tuple[str, str]] = {
    # --- fragment length -------------------------------------------------
    "frag_length": ("frag", "frag_length"),
    "frag_length_bins": ("frag", "frag_length_bins"),
    "frag_length_intervals": ("frag", "frag_length_intervals"),
    # --- coverage --------------------------------------------------------
    "coverage": ("frag", "coverage"),
    "single_coverage": ("frag", "single_coverage"),
    # --- windowed protection score ---------------------------------------
    "wps": ("frag", "wps"),
    "multi_wps": ("frag", "multi_wps"),
    "adjust_wps": ("frag", "adjust_wps"),
    # --- cleavage profile ------------------------------------------------
    "cleavage_profile": ("frag", "cleavage_profile"),
    "multi_cleavage_profile": ("frag", "multi_cleavage_profile"),
    # --- DELFI -----------------------------------------------------------
    "delfi": ("frag", "delfi"),
    "delfi_gc_correct": ("frag", "delfi_gc_correct"),
    "delfi_merge_bins": ("frag", "delfi_merge_bins"),
    # --- end motifs ------------------------------------------------------
    "end_motifs": ("frag", "end_motifs"),
    "region_end_motifs": ("frag", "region_end_motifs"),
    "interval_end_motifs": ("frag", "interval_end_motifs"),
    "EndMotifFreqs": ("frag", "EndMotifFreqs"),
    "EndMotifsIntervals": ("frag", "EndMotifsIntervals"),
    # --- breakpoint motifs -----------------------------------------------
    "breakpoint_motifs": ("frag", "breakpoint_motifs"),
    "region_breakpoint_motifs": ("frag", "region_breakpoint_motifs"),
    "interval_breakpoint_motifs": ("frag", "interval_breakpoint_motifs"),
    "BreakpointMotifFreqs": ("frag", "BreakpointMotifFreqs"),
    "BreakpointMotifsIntervals": ("frag", "BreakpointMotifsIntervals"),
    # --- utilities -------------------------------------------------------
    "frag_generator": ("utils", "frag_generator"),
    "frag_array": ("utils", "frag_array"),
    "frags_in_region": ("utils", "frags_in_region"),
    "frag_bam_to_bed": ("utils", "frag_bam_to_bed"),
    "agg_bw": ("utils", "agg_bw"),
    "filter_file": ("utils", "filter_file"),
    "get_intervals": ("utils", "get_intervals"),
    "overlaps": ("utils", "overlaps"),
    "gen_kmers": ("utils", "gen_kmers"),
    "reverse_complement": ("utils", "reverse_complement"),
    "low_quality_read_pairs": ("utils", "low_quality_read_pairs"),
    "chrom_sizes_to_dict": ("utils", "chrom_sizes_to_dict"),
    "chrom_sizes_to_list": ("utils", "chrom_sizes_to_list"),
    # --- genome annotations ----------------------------------------------
    "GenomeGaps": ("genome", "GenomeGaps"),
    "ContigGaps": ("genome", "ContigGaps"),
    "ucsc_hg19_gap_bed": ("genome", "ucsc_hg19_gap_bed"),
    "b37_gap_bed": ("genome", "b37_gap_bed"),
    "ucsc_hg38_gap_bed": ("genome", "ucsc_hg38_gap_bed"),
    # --- I/O wrappers ----------------------------------------------------
    "ReferenceWrapper": ("io", "ReferenceWrapper"),
    "AlignmentWrapper": ("io", "AlignmentWrapper"),
    "Fragment": ("io", "Fragment"),
}

# Convenience singular aliases requested for the modern namespace.
_ALIASES: dict[str, str] = {
    "end_motif": "end_motifs",
    "breakpoint_motif": "breakpoint_motifs",
}


def __getattr__(name: str):
    """Lazily resolve subpackages and flat-namespace symbols (PEP 562)."""
    import importlib

    if name in _SUBMODULES:
        return importlib.import_module(f".{name}", __name__)

    target = _ALIASES.get(name, name)
    if target in _EXPORTS:
        submodule_name, attr = _EXPORTS[target]
        submodule = importlib.import_module(f".{submodule_name}", __name__)
        value = getattr(submodule, attr)
        globals()[name] = value  # cache for subsequent lookups
        return value

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(
        set(globals())
        | set(_SUBMODULES)
        | set(_EXPORTS)
        | set(_ALIASES)
    )


__all__ = [
    "__version__",
    "exceptions",
    *_SUBMODULES,
    *_EXPORTS,
    *_ALIASES,
    # exception classes are eagerly imported above
    "FinaleToolkitError",
    "InvalidInputError",
    "UnsupportedFormatError",
    "MissingReferenceError",
    "MissingIndexError",
    "ContigNotFoundError",
    "ContigMismatchError",
    "OutOfBoundsError",
]
