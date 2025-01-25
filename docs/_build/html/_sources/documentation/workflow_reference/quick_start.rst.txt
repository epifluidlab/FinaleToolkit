Quick Start
-----------

1.  **Configuration:**  Create a ``params.yaml`` file defining your input, output, and processing options (reference below sections).
2.  **Basic Execution:** Run the workflow with ``snakemake --configfile params.yaml -c <cores> -j <jobs>``.

   * ``-c``: Number of CPU cores to use.
   * ``-j``: Maximum number of concurrent jobs.

3.  **SLURM Execution:** Submit to SLURM to run in the background with ``snakemake --profile slurm_profile > snakemake.log 2>&1 &`` (see ``slurm_profile/config.yaml`` for default settings).