Snakemake Workflow
==================

`finaletoolkit_workflow <https://github.com/epifluidlab/finaletoolkit_workflow>`_
is a `Snakemake <https://snakemake.readthedocs.io>`_ pipeline that automates
cell-free DNA fragmentation feature extraction with FinaleToolkit. It runs any
combination of FinaleToolkit commands over a directory of samples, handles
filtering and per-genome reference setup, and scales from a laptop to a SLURM
cluster.

It supports **hg38** and **T2T-CHM13** out of the box, BED/Fragment/BAM/CRAM
inputs, and parallel processing.

Installation
------------

.. code-block:: bash

   git clone https://github.com/epifluidlab/finaletoolkit_workflow
   cd finaletoolkit_workflow
   conda env create -f environment.yml
   conda activate finaletoolkit_workflow

The environment provides ``finaletoolkit``, ``snakemake``, ``bedtools``,
``htslib``, ``samtools``, and ``pybigwig``.

Genome support
--------------

Two genomes ship with a runnable example config and a one-command reference
setup:

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Genome
     - Example config
     - Setup command
   * - hg38
     - ``params.hg38.yaml``
     - ``scripts/setup_reference.sh hg38 supplement 500``
   * - T2T-CHM13 (hs1)
     - ``params.t2t-chm13.yaml``
     - ``scripts/setup_reference.sh t2t-chm13 supplement 500``

Reference and supplement setup
------------------------------

The workflow needs per-genome supplement files (chrom sizes, ``.2bit``, interval
bins, blacklist, gap, and a mappability bigWig). ``scripts/setup_reference.sh``
builds all of them with the exact filenames the example configs expect:

.. code-block:: bash

   # hg38
   scripts/setup_reference.sh hg38 supplement 500

   # T2T-CHM13
   scripts/setup_reference.sh t2t-chm13 supplement 500

This writes into ``supplement/``::

   <g>.chrom.sizes   <g>.2bit   <g>.<N>kb.bins   <g>.delfi.chrom.sizes
   <g>.blacklist.bed   <g>.45mer.mappability.bw   ( + hg38.gap.bed for hg38 )

where ``g`` is ``hg38`` or ``chm13``. ``<g>.delfi.chrom.sizes`` is restricted to
autosomes plus X and Y, since DELFI requires centromere-bearing contigs.

Running the workflow
--------------------

1. Put input fragment/alignment files in ``input/`` (or set ``input_dir``).
2. Pick a config and run:

.. code-block:: bash

   snakemake --configfile params.t2t-chm13.yaml --cores <N> \
     --rerun-incomplete --default-resources "tmpdir='./tmp'"

``--cores`` sets CPU cores; ``--jobs`` caps concurrent jobs; ``-n`` previews the
DAG without running.

SLURM execution
---------------

The SLURM executor plugin is included. Set your account and partition in
``slurm_profile/config.yaml``, then submit:

.. code-block:: bash

   ./ftk_exc.sh params.t2t-chm13.yaml ./tmp   # runs in the background -> snakemake.log

Mappability filtering
---------------------

Interval bins are kept only if their **mean** mappability over the bin is at
least ``mappability_threshold`` (applied by ``scripts/mappability_filter.py``).
A continuous track is recommended: each position is ``1 / (k-mer occurrences)``,
so the bin mean is the average mappability. The shipped hg38 and T2T-CHM13 tracks
are continuous 45-mer GenMap tracks, fetched automatically by
``setup_reference.sh``.

Configuration
-------------

Configuration is a YAML file (``--configfile``). Every option maps **1:1 to a
FinaleToolkit CLI option** of the same command.

* **Enable a command** with ``<command>: True`` (underscores instead of hyphens,
  e.g. ``adjust-wps`` becomes ``adjust_wps: True``).
* **Set an option** with ``<command>_<option>`` (e.g. ``coverage_mapq: 30``).
* Command **dependencies** are enforced: ``mds`` needs ``end_motifs``,
  ``regional_mds`` needs ``interval_end_motifs``, ``adjust_wps`` needs ``wps``,
  ``agg_bw`` needs ``cleavage_profile``.

Files named in the config (references, blacklists, chrom sizes, bins, refseq,
sites) are resolved inside ``supplement_dir``.

Global settings
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``input_dir``
     - Directory of input samples (default ``input``).
   * - ``output_dir``
     - Output directory (default ``output``).
   * - ``supplement_dir``
     - Directory holding reference/supplement files (default ``supplement``).
   * - ``file_format``
     - ``bed.gz``, ``frag.gz``, ``bam``, or ``cram``.
   * - ``interval_file``
     - BED of intervals/bins used by interval-based commands.
   * - ``mappability_file`` / ``mappability_threshold``
     - Mappability bigWig and minimum mean-mappability cutoff.
   * - ``workers``
     - Global default worker count.

Per-command options
~~~~~~~~~~~~~~~~~~~

Each command's options below carry the FinaleToolkit flag they map to. Optional
``<command>_reference`` is a FASTA needed only for CRAM input. Positional inputs
(``*_refseq_file``, ``*_chrom_sizes``, ``*_bins_file``, ``*_site_bed``,
``*_reference_file``) name supplement files.

**filter-file**: ``filter_file_whitelist`` (``-w``), ``filter_file_blacklist``
(``-b``), ``filter_file_mapq`` (``-q``), ``filter_file_min_length``,
``filter_file_max_length``, ``filter_file_intersect_policy`` (``-p``).

**coverage**: ``coverage_reference`` (``-r``), ``coverage_min_len``
(``--min-length``), ``coverage_max_len`` (``--max-length``),
``coverage_normalize`` (``-n``), ``coverage_scale_factor``,
``coverage_intersect_policy`` (``-p``), ``coverage_mapq`` (``-q``),
``coverage_workers`` (``-t``).

**frag-length-bins**: ``frag_length_bins_reference`` (``-r``),
``frag_length_bins_mapq`` (``-q``), ``frag_length_bins_bin_size``,
``frag_length_bins_policy`` (``-p``), ``frag_length_bins_min_len``,
``frag_length_bins_max_len``, ``frag_length_bins_chrom`` (``-c``),
``frag_length_bins_start`` (``-S``), ``frag_length_bins_end`` (``-E``),
``frag_length_bins_summary_stats``, ``frag_length_bins_short_fraction``
(``--short-threshold``).

**frag-length-intervals**: ``frag_length_intervals_reference`` (``-r``),
``frag_length_intervals_min_len``, ``frag_length_intervals_max_len``,
``frag_length_intervals_policy`` (``-p``), ``frag_length_intervals_short_reads``
(``--short-threshold``), ``frag_length_intervals_mapq`` (``-q``),
``frag_length_intervals_workers`` (``-t``).

**end-motifs** / **interval-end-motifs** / **breakpoint-motifs** /
**interval-breakpoint-motifs**: ``*_refseq_file`` (positional reference),
``*_kmer_length`` (``-k``), ``*_min_len``, ``*_max_len``, ``*_strand``
(``--strand``: ``both``, ``forward``, or ``reverse``), ``*_mapq`` (``-q``),
``*_workers`` (``-t``).

**mds**: ``mds_sep`` (``-s``), ``mds_header`` (``--header``).

**regional-mds**: ``regional_mds_sep`` (``-s``), ``regional_mds_header``
(``--header``), ``regional_mds_miller_madow`` (``--miller-madow``;
Miller-Madow bias correction, off by default).

**wps**: ``wps_site_bed`` (positional REGIONS), ``wps_reference`` (``-r``),
``wps_chrom_sizes`` (``--chrom-sizes``), ``wps_interval_size`` (``-i``),
``wps_window_size`` (``-W``), ``wps_min_len``, ``wps_max_len``, ``wps_mapq``
(``-q``), ``wps_workers`` (``-t``).

**adjust-wps**: ``adjust_wps_chrom_sizes`` (positional), ``adjust_wps_interval_size``
(``-i``), ``adjust_wps_median_window_size`` (``-m``),
``adjust_wps_savgol_window_size``, ``adjust_wps_savgol_poly_deg``,
``adjust_wps_savgol`` (``--savgol``/``--no-savgol``; set ``False`` to disable),
``adjust_wps_mean``, ``adjust_wps_subtract_edges``, ``adjust_wps_edge_size``,
``adjust_wps_workers`` (``-t``).

**delfi**: ``delfi_chrom_sizes``, ``delfi_reference_file``, ``delfi_bins_file``
(positionals), ``delfi_blacklist_file`` (``-b``), ``delfi_gap_file`` or
``delfi_gap_reference_genome`` (``-g``; the latter generates a gap BED via
``gap-bed``), ``delfi_no_gc_correct`` (``--no-gc-correct``), ``delfi_remove_nocov``
(``--remove-nocov``/``--no-remove-nocov``, default ``True``), ``delfi_merge_bins``
(``--merge-bins``/``--no-merge-bins``, default ``True``), ``delfi_merge_size``,
``delfi_mapq`` (``-q``), ``delfi_workers`` (``-t``).

**cleavage-profile**: ``cleavage_profile_chrom_sizes`` (positional),
``cleavage_profile_reference`` (``-r``), ``cleavage_profile_min_len``,
``cleavage_profile_max_len``, ``cleavage_profile_mapq`` (``-q``),
``cleavage_profile_left`` (``--pad-left``), ``cleavage_profile_right``
(``--pad-right``), ``cleavage_profile_workers`` (``-t``).

**agg-bw**: ``agg_bw_median_window_size`` (``-m``), ``agg_bw_mean`` (``--mean``).

.. tip::

   The repository's `params.yaml
   <https://github.com/epifluidlab/finaletoolkit_workflow/blob/main/params.yaml>`_
   is a fully annotated config listing every option with its flag mapping; copy
   it and uncomment what you need.

Output file naming
------------------

* Filtered files get ``.filtered`` before the format
  (e.g. ``sample.filtered.bed.gz``).
* Command outputs insert the command name (e.g. ``sample.frag_length_intervals.bed``).
* Each input is processed for every enabled command.
