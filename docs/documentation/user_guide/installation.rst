==========================
Installation
==========================

-----------------------
Prerequisites
-----------------------
**FinaleToolkit** requires Python ``3.10`` or newer. Some subcommands
(``filter-file``) additionally require ``samtools``, ``bedtools``, ``bgzip``,
and ``tabix`` on your ``PATH``.

----------------------
Installation
----------------------
**FinaleToolkit** is installed with ``pip``. Using a dedicated ``conda``
environment is recommended::

    conda create -n finaletoolkit python=3.11 -y
    conda activate finaletoolkit
    pip install finaletoolkit

To install from a source checkout (for development, editable install)::

    pip install -e .

Verify the installation::

    finaletoolkit --version

------------------
Errors
------------------

If any errors arise during installation, please open a detailed issue in the
GitHub repository.

If you are getting an ``ImportError`` on Mac when running **FinaleToolkit**,
run the following command in the terminal::

    $ brew install curl

If you see a NumPy ``_ARRAY_API not found`` warning, a package in your Python
*user site* (``~/.local``) was compiled against an incompatible NumPy. Run with
``PYTHONNOUSERSITE=1`` or remove the offending package from ``~/.local``.
