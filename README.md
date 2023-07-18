# FinaleTools

Python tools for analysis of cell free DNA fragment data. FinaleTools refers to FragmentatIoN AnaLysis of cEll-free DNA Tools.

This package is currently wip. Here is how to install (very wip):

- Clone repo from https://github.com/epifluidlab/FinaleTools/
- Create a conda environment with all of the dependencies. If only using pip to install dependencies, create venv/conda env and go to next step
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
