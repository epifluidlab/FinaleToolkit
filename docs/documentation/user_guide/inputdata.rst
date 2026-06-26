Input Data
==========

FinaleToolkit works with almost any paired-end sequencing data, in three input
formats. Each one must be indexed so the engine can stream fragments
efficiently.

.. list-table::
   :header-rows: 1
   :widths: 22 30 28 20

   * - Format
     - Extension
     - Required index
     - Reference needed?
   * - BAM
     - ``.bam``
     - ``.bam.bai``
     - No
   * - CRAM
     - ``.cram``
     - ``.cram.crai``
     - Yes (FASTA)
   * - Fragment
     - ``.frag.gz`` or ``.frag.tsv.bgz``
     - ``.tbi`` (Tabix)
     - No

BAM
---

A Binary Alignment Map stores sequence-alignment results: where each read maps
to the reference, plus information about the read itself. It is not
human-readable, but it is space-efficient compared to its plaintext
counterpart, SAM.

.. note::

   BAM input must be BAI-indexed. Keep the ``.bam.bai`` file alongside the
   ``.bam``.

CRAM
----

CRAM is a reference-based binary alignment format. It is smaller than BAM while
carrying the same information.

.. note::

   CRAM input must be CRAI-indexed. Keep the ``.cram.crai`` file alongside the
   ``.cram``.

.. admonition:: CRAM requires a reference FASTA
   :class: important

   Because CRAM encodes reads relative to the reference, you must pass the same
   FASTA reference used during alignment to any subcommand::

      $ finaletoolkit coverage sample.cram intervals.bed -r reference.fa -o cov.bed

   The FASTA needs a ``.fai`` index, which FinaleToolkit creates automatically
   if it is missing. A ``.2bit`` file cannot be used as the reference when
   reading CRAM.

Fragment file
-------------

A Fragment file (``.frag.gz`` or ``.frag.tsv.bgz``) is derived from a BAM file.
It is a block-gzipped BED3+2 file with one row per fragment and these columns:

.. list-table::
   :header-rows: 1
   :widths: 18 18 18 18 18

   * - ``chrom``
     - ``start``
     - ``stop``
     - ``mapq``
     - ``strand``
   * - contig
     - 0-based start
     - end
     - mapping quality
     - ``+`` or ``-``

.. note::

   Fragment files must be Tabix-indexed. Keep the ``.tbi`` file alongside the
   data.

Example::

    #chrom  start   stop    mapq    strand
    chr1    10000   10050   60      +
    chr1    10050   10100   60      -
    chr1    10100   10150   60      +
    chr1    10150   10200   60      -
    chr1    10200   10250   60      +

.. tip::

   Our FinaleDB database hosts ready-to-use fragment files. Learn more at
   `finaledb.research.cchmc.org <http://finaledb.research.cchmc.org>`_.
