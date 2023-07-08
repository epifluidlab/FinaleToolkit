# FinaleTools

Python tools for analysis of cell free DNA fragment data. FinaleTools refers to
FragmentatIoN AnaLysis of cEll-free DNA Tools.

This package is currently wip. Here is how to install (very wip):

- copy from https://github.com/epifluidlab/FinaleTools/
- create a conda environment with all of the dependencies. If only using pip,
    create venv and go to next step
    - for conda, run the following:
    ```
    conda create -n <env_name> -c bioconda -c conda-forge python=3.<version>
    <dependencies>
    ```
    - see pyproject.toml for dependencies
- cd to repo root directory
- run `pip install -e .`
- run `finaletools -h` to see if you did this right