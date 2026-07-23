Features
========

FinaleToolkit computes the established cell-free DNA fragmentomic features
below. Each one is available both as a ``finaletoolkit`` subcommand and as a
Python function. See the :doc:`../cli_reference/index` and
:doc:`../api_reference/index` for exact signatures.

Fragment length
---------------

Fragment-length distributions and per-interval summary statistics across cfDNA
fragments.

:Commands: ``frag-length-bins``, ``frag-length-intervals``

Fragment coverage
-----------------

The number of fragments mapped to a region, reported either raw or normalized
by the sample's total coverage.

:Commands: ``coverage``

Windowed Protection Score (WPS)
-------------------------------

WPS (*Snyder et al. 2016*) quantifies nucleosome protection. For a window
(typically 120 bp) centered on a position, it counts the fragments that fully
span the window and subtracts the fragments that have an endpoint inside it.
The score tracks nucleosome positioning and features such as transcription
start sites and DNase I hypersensitive sites.

:Commands: ``wps``, ``adjust-wps``

DELFI
-----

DELFI (*Cristiano et al. 2019*) detects fragmentation abnormalities in cfDNA.
The score is the ratio of GC-corrected short fragment counts to GC-corrected
long fragment counts. In the original study it was used to categorize cancer
and the associated tissue of origin.

:Commands: ``delfi``

End and breakpoint motifs
-------------------------

Frequencies of k-mer motifs at fragment ends (the cleavage sites) and at
breakpoints. Because these sequences arise from how cfDNA is cut, they can
reflect condition-specific cleavage activity.

:Commands: ``end-motifs``, ``interval-end-motifs``, ``breakpoint-motifs``,
   ``interval-breakpoint-motifs``

Motif Diversity Score (MDS)
---------------------------

MDS (*Jiang et al. 2020*) is the normalized Shannon entropy of the end-motif
k-mer distribution, summarizing end-motif diversity in a single number. The
regional Motif Diversity Score (rMDS) computes that same entropy per genomic
region (*Bandaru et al. 2026*). Because the underlying plug-in entropy
estimator is biased downward in regions with few fragments, ``regional-mds``
accepts an optional ``--miller-madow`` flag that applies a Miller-Madow bias
correction, making rMDS values comparable across regions of differing depth.
It is off by default, so existing output is unchanged unless requested.

:Commands: ``mds`` (whole sample), ``regional-mds`` (per region, optionally
   bias-corrected with ``--miller-madow``)

Cleavage profile
----------------

Cleavage proportion (*Zhou et al. 2022*). At each nucleotide, it is the
fraction of the overlapping fragment ends that fall on that position. The
authors report a relationship to DNA methylation, particularly at CpG sites.

:Commands: ``cleavage-profile``

Intersect policy
----------------

When a feature pulls fragments over an interval ``[start, stop)``, the
``intersect_policy`` argument (CLI: ``-p/--intersect-policy``) decides which
fragments count:

``midpoint`` (default)
   Include a fragment only if its midpoint falls inside the interval. Avoids
   double-counting fragments that straddle a boundary.

``any``
   Include a fragment if it overlaps the interval at all, even by one base pair.

Use ``midpoint`` for binned, genome-wide features; use ``any`` when you care
about every fragment that touches a region.
