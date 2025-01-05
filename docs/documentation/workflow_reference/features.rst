Features
----------

All Finaletoolkit CLI commands are supported *except* for gap-bed and delfi-gc correct. In addition, this workflow supports mappability filtering of interval files, parallelization, slurm, and BED/BAM/CRAM compatability.

Example usage
=============

::

   snakemake --configfile params.yaml -c8 -j 2

This runs the snakemake pipeline using parameters specified in ``params.yaml`` on 8 cores, with up to 2 rules running at a time.

Input files from the specified input directory are processed into the specified output directory. Blacklist and secondary input files for finaletoolkit go in the supplementary directory.

::

   snakemake --profile slurm_profile > snakemake.log 2>&1 &

This runs the snakemake pipeline using flags from ``slurm_profile/config.yaml``. In the sample ``config.yaml`` file given in `this repository <https://github.com/epifluidlab/finaletoolkit_workflow>`_, the command above would create a slurm submission that allows for up to 2GB per snakemake job, up to 4 jobs running in parallel, and up to 8 cores per job. It runs in the background (``&``) with output going to ``snakemake.log``.
