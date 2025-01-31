Quick Start
-----------

*  **Configuration:**  Create a ``params.yaml`` file defining your input, output, and processing options (reference below sections).

*	**Basic Execution:** Once you're in the relevant folder with the ``Snakefile`` present, run the workflow through the following command:

	.. code-block:: bash
 
		cd finaletoolkit_workflow # Enter the folder with the workflow Snakefile

		snakemake --configfile params.yaml --cores <cores> --jobs <jobs>
		# --cores: Number of CPU cores to use.
		# --jobs: Maximum number of concurrent jobs (limited by --cores).

*	**SLURM Execution:** Submit to SLURM to run the workflow through the command below (see ``slurm_profile/config.yaml`` for default settings).

	.. code-block:: bash

		cd finaletoolkit_workflow # Enter the folder with the workflow Snakefile
      
		snakemake  --profile slurm_profile > snakemake.log 2>&1 &
		# Runs the command through params specified in slurm_profile/config.yaml in the background (&),
		# Redirects all command-related output to snakemake.log
