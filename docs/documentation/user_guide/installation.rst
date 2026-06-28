Installation
============

FinaleToolkit needs **Python 3.10 or newer**. The ``filter-file`` command also
needs ``samtools``, ``bedtools``, ``bgzip``, and ``tabix`` on your ``PATH``
(``conda install -c bioconda samtools bedtools htslib``).

Install into a clean environment:

.. tab-set::

   .. tab-item:: conda

      .. code-block:: console

         $ conda create -n finaletoolkit python=3.11 -y
         $ conda activate finaletoolkit
         $ pip install finaletoolkit

   .. tab-item:: venv

      .. code-block:: console

         $ python3 -m venv .venv
         $ source .venv/bin/activate
         $ pip install finaletoolkit

   .. tab-item:: From source

      .. code-block:: console

         $ git clone https://github.com/epifluidlab/finaletoolkit.git
         $ cd finaletoolkit
         $ pip install -e .

Then verify it::

    $ finaletoolkit --version

Troubleshooting
---------------

- **ImportError on macOS:** install a current ``curl`` (``brew install curl``).
- **NumPy ``_ARRAY_API not found``:** a package in ``~/.local`` was built
  against an incompatible NumPy. Run with ``PYTHONNOUSERSITE=1`` or remove it.
- **Anything else:** open an issue on the `tracker
  <https://github.com/epifluidlab/finaletoolkit/issues>`_ with a reproducible
  example, the error output, and your ``finaletoolkit --version``.
