# FinaleTools

Python tools for analysis of cell free DNA fragment data. FinaleTools refers to FragmentatIoN AnaLysis of cEll-free DNA Tools.

## Installation
This package is currently work-in-progress. Here is how to install (wip):

- Clone repo from https://github.com/epifluidlab/FinaleTools/
- Create a conda environment with all of the dependencies. If only using pip to install dependencies, create venv/conda env and go to next step.
  Pip might not install all non-python dependencies.
    - For conda, run the following:
    ```
    conda create -n <env_name> -c bioconda -c conda-forge python=3.<version>
    <dependencies>
    - see pyproject.toml for dependencies
    ```
    - Alternatively if using Linux, you can also use the attached YAML file by navigating to the project directory and running:
      ```
      conda create -f conda_env/environment.yml
      ```
      The resulting environment should be called `samenv3`.
- cd to project directory
- Run `pip install -e .`
- Run `finaletools -h` to see if you did this right

## Quick Start
FinaleTools functions generally accept reads in two file formats:
- Binary Alignment Map (BAM) Files
- FinaleDB Frag Files

Frag files are BED3+2 files with the following format:
`chrom  start  stop  mapq  strand(+/-)`

Frag files can be retrieved from http://finaledb.research.cchmc.org/

Because FinaleTools uses pysam, BAM files should be bai-indexed and Frag files should be tabix-indexed.

To verify FinaleTools has been successfully installed, try
```
$ finaletools -h
usage: finaletools [-h]
                   {coverage,frag-length,frag-length-bins,frag-length-intervals,wps,delfi,filter-bam,adjust-wps,agg-wps,delfi-gc-correct,end-motifs,mds}
                   ...

Calculates fragmentation features given a CRAM/BAM/SAM file

options:
  -h, --help            show this help message and exit

subcommands:
  {coverage,frag-length,frag-length-bins,frag-length-intervals,wps,delfi,filter-bam,adjust-wps,agg-wps,delfi-gc-correct,end-motifs,mds}
```

