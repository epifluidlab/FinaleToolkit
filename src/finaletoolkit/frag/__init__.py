"""
Fragmentation feature extractors.

Re-exports the original public surface unchanged, plus a few additive names
(``multi_cleavage_profile``, ``region_end_motifs``, ``region_breakpoint_motifs``).
"""
from ._frag_length import (
    frag_length,
    frag_length_bins,
    frag_length_intervals,
    FragLengthStats,
)
from ._coverage import coverage, single_coverage, CoverageResult
from ._wps import wps
from ._multi_wps import multi_wps
from ._cleavage_profile import cleavage_profile, multi_cleavage_profile
from ._delfi import delfi
from ._delfi_gc_correct import delfi_gc_correct
from ._delfi_merge_bins import delfi_merge_bins
from ._adjust_wps import adjust_wps
from ._end_motifs import (
    EndMotifFreqs,
    EndMotifsIntervals,
    region_end_motifs,
    end_motifs,
    interval_end_motifs,
)
from ._breakpoint_motifs import (
    BreakpointMotifFreqs,
    BreakpointMotifsIntervals,
    region_breakpoint_motifs,
    breakpoint_motifs,
    interval_breakpoint_motifs,
)

__all__ = [
    "frag_length",
    "frag_length_bins",
    "frag_length_intervals",
    "FragLengthStats",
    "coverage",
    "single_coverage",
    "CoverageResult",
    "wps",
    "multi_wps",
    "cleavage_profile",
    "multi_cleavage_profile",
    "delfi",
    "delfi_gc_correct",
    "delfi_merge_bins",
    "adjust_wps",
    "EndMotifFreqs",
    "EndMotifsIntervals",
    "region_end_motifs",
    "end_motifs",
    "interval_end_motifs",
    "BreakpointMotifFreqs",
    "BreakpointMotifsIntervals",
    "region_breakpoint_motifs",
    "breakpoint_motifs",
    "interval_breakpoint_motifs",
]
