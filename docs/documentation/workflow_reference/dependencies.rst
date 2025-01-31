Dependencies
------------

This workflow relies on the following tools being installed and accessible by your system PATH. We recommend that you install FinaleToolkit through ``pip`` and the other packages through ``conda`` in the Bioconda channel:

* ``finaletoolkit``: A command-line tool for epigenomic feature extraction.
* ``snakemake``: A workflow engine that determines which operations ("rules") to carry out on genomic files.
* ``bedtools``: A suite of utilities for working with and manipulating genomic intervals.
* ``htslib``: A library that includes ``bgzip```, necessary to GZIP uncompressed BED files. 
* ``samtools``: A set of tools for manipulating and analyzing sequencing BAM/CRAM data