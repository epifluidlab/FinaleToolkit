CLI Reference
=============

.. tip::

   Every subcommand also prints its own options and a worked example with
   ``finaletoolkit <command> --help``. Flags are consistent across commands:
   ``-o/--output``, ``-r/--reference``, ``-q/--min-mapq``,
   ``--min-length`` / ``--max-length``, ``-t/--threads``, ``-v/--verbose``,
   and ``-k/--kmer-length``.

.. click:: finaletoolkit.cli.main_cli:main_cli
   :prog: finaletoolkit
   :nested: full
