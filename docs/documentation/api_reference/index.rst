API Reference
=============

Every public feature is reachable directly from the top-level ``finaletoolkit``
namespace (see :doc:`namespace`). You can also import each feature from its
submodule (``finaletoolkit.frag``, ``finaletoolkit.utils``,
``finaletoolkit.genome``, ``finaletoolkit.io``).

.. code-block:: python

   import finaletoolkit as ftk

   cov = ftk.coverage("sample.bam", "intervals.bed", output_file=None)
   gaps = ftk.GenomeGaps.hg38()

.. rubric:: Contents

* :doc:`namespace` lists everything importable from ``finaletoolkit``.
* :doc:`basicfeatures` covers coverage and fragment length.
* :doc:`wps` covers the Windowed Protection Score.
* :doc:`delfi` covers DELFI ratios and GC correction.
* :doc:`endmotifs` covers end-motif frequencies and MDS.
* :doc:`breakpointmotifs` covers breakpoint-motif frequencies.
* :doc:`cleavageprofile` covers the per-nucleotide cleavage proportion.
* :doc:`fragfile` covers filtering, aggregation, and helpers.
* :doc:`genomeutils` covers genome gaps and annotation tracks.
* :doc:`io` covers reference and alignment access.
* :doc:`exceptions` covers the typed exception hierarchy.

.. toctree::
   :hidden:
   :maxdepth: 1

   namespace
   basicfeatures
   wps
   delfi
   endmotifs
   breakpointmotifs
   cleavageprofile
   fragfile
   genomeutils
   io
   exceptions
