===
API
===

Basic Features
==============

.. autofunction:: finaletools.frag.coverage

.. autofunction:: finaletools.frag.frag_length

.. autofunction:: finaletools.frag.frag_length_bins

Window Protection Score
=======================

.. autofunction:: finaletools.frag.wps

.. autofunction:: finaletools.frag.multi_wps

.. autofunction:: finaletools.frag.adjust_wps

.. autofunction:: finaletools.frag.agg_wps

DELFI
=====

.. autofunction:: finaletools.frag.delfi

.. autofunction:: finaletools.frag.delfi_gc_correct

.. autofunction:: finaletools.frag.delfi_merge_bins

End-motifs
==========

.. autoclass:: finaletools.frag.EndMotifFreqs
   :members:

.. autoclass:: finaletools.frag.EndMotifsIntervals
   :members:

.. autofunction:: finaletools.frag.region_end_motifs

.. autofunction:: finaletools.frag.end_motifs

.. autofunction:: finaletools.frag.interval_end_motifs

Cleavage Profile
================
.. autofunction:: finaletools.frag.cleavage_profile


Frag File Utilities
===================

.. autofunction:: finaletools.utils.filter_bam

.. autofunction:: finaletools.utils.genome2list

Genome Utilities
================
.. autoclass:: finaletools.genome.GenomeGaps
   :members:

.. autoclass:: finaletools.genome.ContigGaps
   :members:

.. autofunction:: finaletools.genome.ucsc_hg19_gap_bed

.. autofunction:: finaletools.genome.b37_gap_bed

.. autofunction:: finaletools.genome.ucsc_hg38_gap_bed