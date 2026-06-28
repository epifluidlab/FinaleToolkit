"""
Utility layer: fragment I/O helpers, interval math, k-mer helpers, file
filtering, bigWig aggregation, and validation.

The public surface is identical to the original package; validation helpers are
re-exported additively.
"""
from .utils import (
    chrom_sizes_to_list,
    chrom_sizes_to_dict,
    frag_bam_to_bed,
    frags_in_region,
    frag_array,
    low_quality_read_pairs,
    _not_read1_or_low_quality,
    get_intervals,
    overlaps,
    gen_kmers,
    reverse_complement,
)
from ._intervals import (
    _merge_overlapping_intervals,
    _reduce_overlaps_in_file,
    _convert_to_list,
    _merge_all_intervals,
)
from ._comparison import _none_eq, _none_geq, _none_leq
from ._agg_bw import agg_bw
from ._filter_file import filter_file
from ._frag_generator import frag_generator
from .validation import validate_compatible_contigs, valid_interval

__all__ = [
    "chrom_sizes_to_list",
    "chrom_sizes_to_dict",
    "frag_bam_to_bed",
    "frags_in_region",
    "frag_array",
    "low_quality_read_pairs",
    "_not_read1_or_low_quality",
    "get_intervals",
    "overlaps",
    "gen_kmers",
    "reverse_complement",
    "_merge_overlapping_intervals",
    "_reduce_overlaps_in_file",
    "_convert_to_list",
    "_merge_all_intervals",
    "_none_eq",
    "_none_geq",
    "_none_leq",
    "agg_bw",
    "filter_file",
    "frag_generator",
    "validate_compatible_contigs",
    "valid_interval",
]
