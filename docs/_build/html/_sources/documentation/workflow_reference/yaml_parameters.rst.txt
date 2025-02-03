YAML Parameters
---------------

**Required:**
    *   ``input_dir``: Path to the input directory. Defaults to ``input`` if not specified.
    *   ``output_dir``: Path to the output directory. Defaults to ``output`` if not specified.
    *   ``file_format``: ``"bed.gz"``, ``"frag.gz"``, ``"bam"``, or ``"cram"`` indicating the format of the input files. Defaults to ``bed.gz`` if not specified.

**Optional:**
    *   ``supplement_dir``: Path to supplemental files directory. Defaults to ``supplement`` if not specified.
    *   ``mappability_file``: Name of the bigWig mappability file in ``supplement_dir``.
    *    ``mappability_threshold``: Minimum average mappability score (0.0-1.0) for interval filtering.
    *   ``interval_file``: Path to interval file in ``supplement_dir``.
    *   ``finaletoolkit_command: True/False``: Enables a specific Finaletoolkit command, using hyphens replaced by underscores (e.g., ``adjust-wps`` becomes ``adjust_wps: True``).
    *   ``finaletoolkit_command_flag: value``: Sets flags for a Finaletoolkit command (e.g., ``adjust_wps_max_length: 250``). Flags that take input files, output files, or ``verbose`` flags do not exist here.  ``mapping_quality`` is shortened to ``mapq`` for flags (e.g., ``coverage_mapping_quality`` becomes ``coverage_mapq``).