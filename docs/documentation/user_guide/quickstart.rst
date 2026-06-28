Quickstart
==========

FinaleToolkit works the same way from the command line or from Python. Pick
whichever fits your workflow. The features and results are identical.

Command line
------------

List the available subcommands::

    $ finaletoolkit --help

A few representative runs:

.. code-block:: console

   $ finaletoolkit coverage sample.bam intervals.bed -o coverage.bed
   $ finaletoolkit wps sample.bam tss.bed --chrom-sizes hg38.chrom.sizes -o wps.bw -t 8
   $ finaletoolkit end-motifs sample.bam hg38.2bit -k 4 -o motifs.tsv
   $ finaletoolkit delfi sample.bam autosomes.chrom.sizes hg19.2bit bins.bed -g hg19 -o delfi.tsv

Every subcommand has its own ``--help`` with a worked example. Flags are
consistent across commands, so once you learn them they transfer everywhere:

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Flag
     - Meaning
   * - ``-o`` / ``--output``
     - Output file path. Use ``-`` to write to standard output (stdout).
   * - ``-r`` / ``--reference``
     - Reference FASTA file (required for CRAM input).
   * - ``-q`` / ``--min-mapq``
     - Minimum mapping quality.
   * - ``--min-length`` / ``--max-length``
     - Fragment-length bounds, in base pairs.
   * - ``-t`` / ``--threads``
     - Number of worker processes.
   * - ``-v`` / ``--verbose``
     - Increase verbosity. Repeat for more detail (``-vv``).
   * - ``-k`` / ``--kmer-length``
     - k-mer length (motif commands).

.. tip::

   ``-`` as an output means "write to standard output" instead of a file, so
   you can pipe results into another tool, for example
   ``finaletoolkit mds motifs.tsv -o - | less``.

See the :doc:`../cli_reference/index` for the complete reference.

Python API
----------

Everything is reachable from the top-level ``finaletoolkit`` namespace:

.. code-block:: python

   import finaletoolkit as ftk

   cov = ftk.coverage("sample.bam", "intervals.bed", output_file=None)
   motifs = ftk.end_motifs("sample.bam", "hg38.2bit", k=4)
   mds = motifs.motif_diversity_score()

The package is also organized into submodules, which remain importable:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Submodule
     - Contents
   * - ``frag``
     - Fragmentomic feature generation.
   * - ``genome``
     - Utilities for genome tracks and gaps.
   * - ``utils``
     - Helpers that simplify feature generation.
   * - ``io``
     - Reference and alignment or fragment wrappers.
   * - ``cli``
     - Command-line interface.

To load a specific function from its submodule::

    >>> from finaletoolkit.frag import delfi

.. seealso::

   For end-to-end tutorials, see the
   `wiki <https://github.com/epifluidlab/FinaleToolkit/wiki>`_. For the full
   Python surface, see the :doc:`../api_reference/index`.
