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

.. _delfi-bins:

**Bins.** To replicate the original study, use the 100kb hg19 HiC A/B
compartment bins the DELFI scripts load:
`hic_compartments_100kb_ebv_2014.txt
<https://github.com/Jfortin1/HiC_AB_Compartments/blob/master/data/hic_compartments_100kb_ebv_2014.txt>`_
from `HiC_AB_Compartments <https://github.com/Jfortin1/HiC_AB_Compartments>`_
(Fortin & Hansen 2015). That file is space-delimited with an uncommented
header, so convert it to a tab-delimited BED with a ``#`` header first:

.. code-block:: console

   $ curl -sL https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt \
       | awk 'BEGIN{OFS="\t"} NR==1{$1="#"$1} {$1=$1; print}' > hic_100kb.bed

Any BED of uniform bins also works. For example, 100kb autosomal bins for
hg38 (requires ``bedtools``), which also writes the autosome-only
``CHROM_SIZES`` file to pass to ``delfi``:

.. code-block:: console

   $ curl -s https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes \
       | awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' | sort -k1.4,1n > hg38.autosomes.chrom.sizes
   $ bedtools makewindows -g hg38.autosomes.chrom.sizes -w 100000 \
       | awk '$3 - $2 == 100000' > hg38_100kb_bins.bed
   $ finaletoolkit delfi sample.bam hg38.autosomes.chrom.sizes hg38.2bit hg38_100kb_bins.bed \
       -g hg38 --no-remove-nocov -o delfi.tsv

The final ``awk`` drops the shorter last bin of each chromosome (``delfi``
would otherwise drop it with a warning). When bins are merged (the default), every
bin must be the same size and ``--merge-size`` (default 5Mb) must be a multiple
of it; a shorter last bin on a contig, as ``bedtools makewindows`` produces, is
dropped with a warning. Bins are merged by count (``merge size / bin size``
bins per window) within each chromosome arm, so a window can span more than
``--merge-size`` where bins were removed. Without a gap file (``-g``), arms are
unknown and each contig is merged as a whole.

**Coordinates.** Both half-open BED bins (``stop - start = 100000``) and
end-inclusive bins (``stop - start = 99999``, as in the HiC file above) are
accepted, and the bin size is inferred accordingly. FinaleToolkit reads bin
coordinates as written (0-based, half-open). The original R scripts read the
HiC file as 1-based and closed, so each bin there covers 1 bp earlier and
1 bp more than it does here (``chr1 700000 799999`` is ``[699999, 799999)``
in R and ``[700000, 799999)`` in FinaleToolkit). The effect on fragment
counts is negligible, but results will not match the R output to the
fragment.

**No-coverage bins.** ``--remove-nocov`` (on by default) drops the two bins the
original scripts found to have no coverage. They are identified by **row
position** (the 8,780th and 13,665th bins left after gap and arm filtering),
exactly as in the original scripts, so the option is only meaningful with the
HiC bins above, the hg19 gap file (``-g hg19``), and all autosomes in
``CHROM_SIZES``. With any other bins file, reference, or subset of
chromosomes it removes two arbitrary bins; pass ``--no-remove-nocov``.

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
