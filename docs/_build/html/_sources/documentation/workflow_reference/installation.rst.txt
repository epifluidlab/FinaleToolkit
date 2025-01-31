Installation
------------

Reference the following command for setup and execution:

.. code-block:: console

    $ git clone https://github.com/epifluidlab/finaletoolkit_workflow # Download the repository containing the workflow
    $ cd finaletoolkit_workflow # Enter the repository folder
    $ conda env create -f environment.yml # Create environment with relevant conda packages
    $ conda activate finaletoolkit_workflow # Use environment for finaletoolkit-workflow
    $ pip install finaletoolkit # Install finaletoolkit seperately through pip
    $ snakemake --configfile params.yaml --cores 4 --jobs 2 # Run with parameters set in params.yaml