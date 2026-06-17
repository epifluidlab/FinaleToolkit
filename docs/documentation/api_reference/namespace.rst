Top-Level Namespace
======================================

Everything below is importable directly from ``finaletoolkit`` without knowing
the internal submodule layout::

    import finaletoolkit as ftk

    results = ftk.coverage("sample.bam", "intervals.bed", output_file=None)
    motifs = ftk.end_motifs("sample.bam", "hg38.2bit", k=4)
    gaps = ftk.GenomeGaps.hg38()

The legacy import paths remain valid, so existing code keeps working::

    from finaletoolkit.frag import wps, multi_wps, delfi
    from finaletoolkit.utils import frag_generator, agg_bw, filter_file
    from finaletoolkit.genome import GenomeGaps, ContigGaps
    from finaletoolkit.io import ReferenceWrapper, AlignmentWrapper, Fragment

Fragmentation features
----------------------

- ``frag_length``, ``frag_length_bins``, ``frag_length_intervals``
- ``coverage``, ``single_coverage``
- ``wps``, ``multi_wps``, ``adjust_wps``
- ``cleavage_profile``, ``multi_cleavage_profile``
- ``delfi``, ``delfi_gc_correct``, ``delfi_merge_bins``
- ``end_motifs``, ``region_end_motifs``, ``interval_end_motifs``
- ``breakpoint_motifs``, ``region_breakpoint_motifs``,
  ``interval_breakpoint_motifs``
- ``end_motif`` / ``breakpoint_motif`` (aliases of the plural functions)

Result and container types
--------------------------

- ``CoverageResult``, ``FragLengthStats``
- ``EndMotifFreqs``, ``EndMotifsIntervals``
- ``BreakpointMotifFreqs``, ``BreakpointMotifsIntervals``

Utilities
---------

- ``frag_generator``, ``frag_array``, ``frags_in_region``, ``frag_bam_to_bed``
- ``agg_bw``, ``filter_file``
- ``get_intervals``, ``overlaps``, ``gen_kmers``, ``reverse_complement``
- ``chrom_sizes_to_dict``, ``chrom_sizes_to_list``, ``low_quality_read_pairs``

I/O and genome annotations
--------------------------

- ``ReferenceWrapper``, ``AlignmentWrapper``, ``Fragment``
- ``GenomeGaps``, ``ContigGaps``
- ``ucsc_hg19_gap_bed``, ``b37_gap_bed``, ``ucsc_hg38_gap_bed``

Exceptions
----------

See :doc:`exceptions`.
