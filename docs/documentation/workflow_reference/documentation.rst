Documentation
-------------

Basic I/O 
=========

The ``Snakefile`` in `this repository <https://github.com/epifluidlab/finaletoolkit_workflow>`_ and a YAML file specifying the parameters to process your data must be present when running the workflow, in a manner similar to the example usage section above.

YAML parameters (mappings) are of the form: ``key: value``, where values that are strings must be wrapped in double quotes, and each parameter is on its own line. This workflow requires the keys ``input_dir`` and ``output_dir`` to specify the folders in the working directory of which genomic data should be taken from and where it should be processed into.

``supplement_dir`` is also necessary to denote the folder of files that are not processed themselves in a command, but are needed to process input data. For example, you would put your blacklist, whitelist, mappability, and interval files in the folder specified by ``supplement_dir``.

Files in the ``output_dir`` directory will have an extension added before their format type. All files that make their way from the input to the output directory (after being run through ``filter-file`` if present) will have a ``.final`` extension before the type denoting the format. For example, ``input/file.bed.gz`` would be filtered into ``output/file.final.bed.gz``. 

Additionally if Finaletoolkit commands other than ``filter-file`` are specified, the secondary extension will be the name of the command the file was run through with hyphens substituted for underscores. For example, if you run ``frag-length-bins``, files in an input directory like ``input/file.bed.gz`` would be run through Finaletoolkit and outputted as ``output/file.frag_length_bins.bed.gz``. If multiple Finaletoolkit commands are specified, then each file will be outputted for each Finaletoolkit command.

Filtering by mappability
========================

``mappability_file`` and ``mappability_threshold`` are used to filter interval files located in ``supplement_dir``. ``mappability_file`` should be the name of the bigWig file in ``supplement_dir`` that filters all interval files, and ``mappability_threshold`` should be a floating point value from 0.0 to 1.0 to denote the minimum average mappability value that intervals will be filtered by as per the mappability file. Finaletoolkit commands that take in these interval files will recieve the interval file filtered by mappability.

Using Finaletoolkit in this workflow
====================================

Finaletoolkit CLI commands and flags directly correspond to parameters in this workflow, with hyphens turning into underscores, and flags being seperated from their command by an underscore. For example, ``adjust-wps`` would be ``adjust_wps``, and the max length flag of this command would be set through the key ``adjust_wps_max_length``.

If you want to use a Finaletoolkit command, set the converted command name in YAML to ``True``. For example, if you wanted to run ``adjust-wps`` on your input files, include the line ``adjust_wps: True``. Flags may be set through the naming scheme as specified in the paragraph above (``adjust_wps_max_length: 250``). Flags for this workflow do not exist for ``input_file``, ``output_file``, deperecated flags, or ``verbose``, which is on by default.

Certain commands require prerequisite commands to run. For example, you must run ``wps`` to run ``adjust-wps`` since the ``adjust-wps`` rule takes in the output of ``wps``. Similarly, ``mds`` requires ``end-motifs``, ``interval-end-motifs`` requires ``interval-mds``, and ``agg-bw`` requires ``cleavage-profile``.

Filter-file
===========

The ``filter-file`` command can be used in a different manner than the other Finaletoolkit commands. If only ``filter-file`` is set to run, then it will be the output of the workflow. However, if other Finaletoolkit commands are set to run, then they will take in the output of ``filter-file``.

Files associated with ``filter_file_blacklist_file`` and ``filter_file_whitelist_file`` will automatically be block gzipped and indexed if they do not end in ``.gz``.