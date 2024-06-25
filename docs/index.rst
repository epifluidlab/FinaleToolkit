
FinaleToolkit
=========================================

-------
About
-------

FinaleToolkit (FragmentatIoN AnaLysis of cEll-free DNA Toolkit) is a package and standalone program to extract fragmentation features of cell-free DNA from paired-end sequencing data.

- **Website:** https://epifluidlab.github.io/finaletoolkit-docs/
- **PyPI:** https://pypi.org/project/finaletoolkit/
- **Documentation:** https://epifluidlab.github.io/finaletoolkit-docs/documentation
- **Source code:** https://github.com/epifluidlab/finaletoolkit
- **Bug reports:** https://github.com/epifluidlab/finaletoolkit/issues


------------------
Documentation
------------------
.. toctree::
   :maxdepth: 2
   
   documentation/index
   documentation/user_guide/index
   documentation/cli_reference/index
   documentation/api_reference/index
   
-----------------
Citations
-----------------

If FinaleToolkit is integral to a scientific publication, please cite it. A paper describing FinaleToolkit has been written. Here is a ready-made BibTeX entry::

    @article {Li2024.05.29.596414,
	author = {Li, James W. and Bandaru, Ravi and Liu, Yaping},
	title = {FinaleToolkit: Accelerating Cell-Free DNA Fragmentation Analysis with a High-Speed Computational Toolkit},
	elocation-id = {2024.05.29.596414},
	year = {2024},
	doi = {10.1101/2024.05.29.596414},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Cell-free DNA (cfDNA) fragmentation pattern represents a promising non-invasive biomarker for disease diagnosis and prognosis. Numerous fragmentation features, such as end motif and window protection score (WPS), have been characterized in cfDNA genomic sequencing. However, the analytical tools developed in these studies are often not released to the liquid biopsy community or are inefficient for genome-wide analysis in large datasets. To address this gap, we have developed FinaleToolkit, a fast and memory efficient Python package designed to generate comprehensive fragmentation features from large cfDNA genomic sequencing data. For instance, FinaleToolkit can generate genome-wide WPS features from a \~{}100X cfDNA whole-genome sequencing (WGS) dataset in 1.2 hours using 16 CPU cores, offering up to a \~{}50-fold increase in processing speed compared to original implementations in the same dataset. We have benchmarked FinaleToolkit against original studies or implementations where possible, confirming its efficacy. Furthermore, FinaleToolkit enabled the genome-wide analysis of fragmentation patterns over arbitrary genomic intervals, significantly boosting the performance for cancer early detection. FinaleToolkit is open source and thoroughly documented with both command line interface and Python application programming interface (API) to facilitate its widespread adoption and use within the research community: https://github.com/epifluidlab/FinaleToolkitCompeting Interest StatementY.L. owns stocks from Freenome Inc. The remaining authors declare no competing interests.},
	URL = {https://www.biorxiv.org/content/early/2024/06/02/2024.05.29.596414},
	eprint = {https://www.biorxiv.org/content/early/2024/06/02/2024.05.29.596414.full.pdf},
	journal = {bioRxiv}
    }


-----------------
Contact
-----------------

| **James Wenhan Li**
| Email: lijw21@wfu.edu

| **Ravi Bandaru**
| Email: ravi.bandaru@northwestern.edu

| **Yaping Liu**
| Email: yaping@northwestern.edu

-----------------
License
-----------------

Please refer to the `MIT license <https://github.com/epifluidlab/FinaleToolkit/blob/main/LICENSE>`_.