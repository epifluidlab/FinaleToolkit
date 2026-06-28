.. raw:: html

   <div class="ftk-hero">
     <h1>FinaleToolkit</h1>
     <p class="ftk-tagline">
       Extract cell-free DNA fragmentation features from whole-genome
       sequencing, through one command line and Python API.
     </p>
   </div>

**FinaleToolkit** (*FragmentatIoN AnaLysis of cEll-free DNA Toolkit*) computes
established cfDNA fragmentomic features (coverage, fragment length, WPS, DELFI,
motifs, MDS, cleavage profiles) on a fast, streaming engine.

.. tab-set::

   .. tab-item:: Command line

      .. code-block:: console

         $ finaletoolkit coverage sample.bam intervals.bed -o coverage.bed
         $ finaletoolkit wps sample.bam tss.bed --chrom-sizes hg38.chrom.sizes -o wps.bw -t 8
         $ finaletoolkit end-motifs sample.bam hg38.2bit -k 4 -o motifs.tsv

   .. tab-item:: Python

      .. code-block:: python

         import finaletoolkit as ftk

         cov = ftk.coverage("sample.bam", "intervals.bed", output_file=None)
         motifs = ftk.end_motifs("sample.bam", "hg38.2bit", k=4)
         mds = motifs.motif_diversity_score()

Features
--------

.. raw:: html

   <div class="ftk-features">
     <div class="ftk-feature">
       <div class="ftk-feature-name">Fragment length</div>
       <div class="ftk-cmds"><span class="ftk-cmd">frag-length-bins</span><span class="ftk-cmd">frag-length-intervals</span></div>
     </div>
     <div class="ftk-feature">
       <div class="ftk-feature-name">Coverage</div>
       <div class="ftk-cmds"><span class="ftk-cmd">coverage</span></div>
     </div>
     <div class="ftk-feature">
       <div class="ftk-feature-name">Windowed Protection Score</div>
       <div class="ftk-cmds"><span class="ftk-cmd">wps</span><span class="ftk-cmd">adjust-wps</span></div>
     </div>
     <div class="ftk-feature">
       <div class="ftk-feature-name">DELFI</div>
       <div class="ftk-cmds"><span class="ftk-cmd">delfi</span></div>
     </div>
     <div class="ftk-feature">
       <div class="ftk-feature-name">End &amp; breakpoint motifs</div>
       <div class="ftk-cmds"><span class="ftk-cmd">end-motifs</span><span class="ftk-cmd">breakpoint-motifs</span></div>
     </div>
     <div class="ftk-feature">
       <div class="ftk-feature-name">Motif Diversity Score</div>
       <div class="ftk-cmds"><span class="ftk-cmd">mds</span><span class="ftk-cmd">regional-mds</span></div>
     </div>
     <div class="ftk-feature">
       <div class="ftk-feature-name">Cleavage profile</div>
       <div class="ftk-cmds"><span class="ftk-cmd">cleavage-profile</span></div>
     </div>
   </div>

Links
-----

.. raw:: html

   <div class="ftk-links">
     <a href="https://github.com/epifluidlab/finaletoolkit">GitHub</a>
     <a href="https://github.com/epifluidlab/finaletoolkit_workflow">Workflow</a>
     <a href="https://pypi.org/project/finaletoolkit/">PyPI</a>
     <a href="http://finaledb.research.cchmc.org">FinaleDB</a>
     <a href="https://github.com/epifluidlab/FinaleToolkit/wiki">Wiki</a>
     <a href="https://github.com/epifluidlab/finaletoolkit/issues">Issues</a>
   </div>

Citation
--------

.. raw:: html

   <div class="ftk-cite">
     <div class="ftk-cite-title">FinaleToolkit: Accelerating Cell-Free DNA Fragmentation Analysis with a High-Speed Computational Toolkit</div>
     <div class="ftk-cite-meta">Li, Bandaru, Baliga &amp; Liu &middot; Bioinformatics Advances, 2025</div>
     <div class="ftk-cite-actions">
       <a href="https://doi.org/10.1093/bioadv/vbaf236">doi:10.1093/bioadv/vbaf236</a>
     </div>
   </div>

.. toctree::
   :hidden:
   :caption: User Guide

   documentation/user_guide/installation
   documentation/user_guide/quickstart
   documentation/user_guide/inputdata
   documentation/user_guide/features

.. toctree::
   :hidden:
   :caption: Reference

   documentation/cli_reference/index
   documentation/api_reference/index

.. toctree::
   :hidden:
   :caption: Workflow

   documentation/workflow/index
