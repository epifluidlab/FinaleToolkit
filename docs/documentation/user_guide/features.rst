
Features
=========================================

**FinaleToolkit** has support for the following cell-free DNA features:


-----------------------
Fragment Length
-----------------------

**FinaleToolkit** can calculate the fragment length distribution and summary statistics of cell-free DNA fragments.

-----------------------
Fragment Coverage
-----------------------

**FinaleToolkit** can calculate fragment coverage (the number of reads mapped to a specific region of the genome). This coverage can be normalized by the total number of reads in the sample, or simply the raw value.

---------------------------------
Windowed Protection Score (WPS)
---------------------------------

**FinaleToolkit** can calculate the Windowed Protection Score (WPS), which is a metric designed by *Snyder et al. 2016*.

It is used to quantify the level of protection or coverage of DNA fragments across a specific genomic region. It's calculated by assessing the number of DNA fragments that fully span a defined window (typically 120 base pairs) centered around a particular genomic coordinate, and then subtracting the number of fragments with an endpoint falling within the same window. This score is correlated with the location of nucleosomes and other genomic features, including transcriptional start sites (TSS) and DNase I hypersensitive sites (DHSs).

-----------------------
DELFI
-----------------------

**FinaleToolkit** can calculate DELFI, which is a metric designed by *Cristiano et al. 2019*. 

DELFI is a metric that was introduced to identify abnormalities in cfDNA from their fragmentation patterns. In the original paper, DELFI was used to categorize patients with cancer and the associated tumor tissues of origin. The associated DELFI score is the ratio between the GC%-corrected short fragment count and the GC%-corrected long fragment count.

-----------------------
End Motifs
-----------------------

**FinaleToolkit** can calculate the frequency of end-motif k-mers of cell-free DNA fragments. Since end motifs are specific sequences found at the ends of cfDNA fragments resulting from cleavage, they can be used to potentially detect patterns associated with certain conditions or diseases.

-----------------------------
Motif Diversity Score (MDS)
-----------------------------

**FinaleToolkit** can calculate the Motif Diversity Score (MDS), which is a metric designed by *Jiang et al. 2020*.

The MDS is a metric that quantifies the diversity of cfDNA end motifs. It is the normalized Shannon entropy of the categorical distribution of all possible end-motif k-mers.

-----------------------
Cleavage Profile
-----------------------

**FinaleToolkit** can calculate the Cleavage Proportion, which is a metric designed by *Zhou et al. 2022*.

The cleavage proportion is, by definition, the number of fragment ends that fall within a nucleotide location over the total number of fragment ends that overlap that nucleotide location. According to the authors, this metric shows a relationship to DNA methylation, specifically methylation at CpG sites.