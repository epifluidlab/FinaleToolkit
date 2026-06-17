
Quickstart
=========================================

------------------------
Command line
------------------------

**FinaleToolkit** is intended to be run directly from a terminal. List the
available subcommands with::

    $ finaletoolkit --help

A few examples::

    $ finaletoolkit coverage sample.bam intervals.bed -o coverage.bed
    $ finaletoolkit wps sample.bam tss.bed --chrom-sizes hg38.chrom.sizes -o wps.bw -t 8
    $ finaletoolkit end-motifs sample.bam hg38.2bit -k 4 -o motifs.tsv
    $ finaletoolkit delfi sample.bam autosomes.chrom.sizes hg19.2bit bins.bed -g hg19 -o delfi.tsv

Every subcommand has ``--help`` with an example invocation. The flag scheme is
consistent across subcommands: ``-o/--output``, ``-r/--reference``,
``-q/--min-mapq``, ``--min-length`` / ``--max-length``, ``-t/--threads``,
``-v/--verbose``, ``-k/--kmer-length``. See the :doc:`../cli_reference/index`
for the full reference.

------------------------
Python API
------------------------

Everything is reachable from the top-level ``finaletoolkit`` namespace::

    >>> import finaletoolkit as ftk
    >>> cov = ftk.coverage("sample.bam", "intervals.bed", output_file=None)
    >>> motifs = ftk.end_motifs("sample.bam", "hg38.2bit", k=4)
    >>> motifs.motif_diversity_score()

The package is also organized into submodules, which remain importable:

* ``frag`` -- fragmentomic feature generation
* ``genome`` -- utilities for genome tracks
* ``utils`` -- helpers that simplify feature generation
* ``io`` -- reference and alignment wrappers
* ``cli`` -- command line interface

To load a specific function from its submodule::

    >>> from finaletoolkit.frag import delfi

For more detailed tutorials, please check out our
`wiki <https://github.com/epifluidlab/FinaleToolkit/wiki>`_.
