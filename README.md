# <img alt="dna with letters FT" src="https://github.com/epifluidlab/FinaleToolkit/blob/b99b38e22a3b07ee9b7e0fa44a488a1eb5442efe/docs/_static/finaletoolkit_logo_rounded.png?raw=true" height="60"> ‎ ‎ ‎FinaleToolkit
<summary><h3>Table of Contents</h2></summary>
<ol>
  <li><a href="#about-the-project">About The Project</a></li>
  <li><a href="#installation">Installation</a></li>
  <li>
    <a href="#usage">Usage</a>
    <ul>
      <li><a href="#functionality">Functionality</a></li>
      <li><a href="#documentation">Documentation</a></li>
      <li><a href="#wikitutorials">Wiki/Tutorials</a></li>
      <li><a href="#compatible-file-formats">Compatible File Formats</a></li>
      <li><a href="#using-fragment-files">Using Fragment Files</a></li>
    </ul>
  </li>
  <li><a href="#contact">Contact</a></li>
  <li><a href="#license">License</a></li>
</ol>




## About The Project
FinaleToolkit (**F**ragmentat**I**o**N** **A**na**L**ysis of c**E**ll-free DNA 
**Toolkit**) is a package and standalone program to extract fragmentation
features of cell-free DNA from paired-end sequencing data.

### Citation
If you use FinaleToolkit in your research, please consider citing our paper:

Li J*, Bandaru R*, Liu Y (2024) FinaleToolkit: Accelerating Cell-Free DNA Fragmentation Analysis with a High-Speed Computational Toolkit. BioRxiv Preprint [![Static Badge](https://img.shields.io/badge/DOI-10.1101%2F2024.05.29.596414-blue?style=flat-square)](https://doi.org/10.1101/2024.05.29.596414)


## Installation
You can install the package using `pip`.
```
$ pip install finaletoolkit
```

## Usage

### Functionality

FinaleToolkit has support for the following cell-free DNA fragmentation features:

- Fragment Length
- Coverage
- End Motifs [![DOI](https://img.shields.io/badge/DOI-10.1158%2F2159--8290.CD--19--0622-blue?style=flat-square&label=DOI)](https://doi.org/10.1158/2159-8290.cd-19-0622)
- Motif Diversity Score [![DOI](https://img.shields.io/badge/DOI-10.1158%2F2159--8290.CD--19--0622-blue?style=flat-square)](https://doi.org/10.1158/2159-8290.CD-19-0622)
- Windowed Protection Score [![DOI](https://img.shields.io/badge/DOI-110.1016%2Fj.cell.2015.11.050-blue?style=flat-square)](https://doi.org/10.1016/j.cell.2015.11.050)
- DELFI [![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41586--019--1272--6-blue?style=flat-square&link=https%3A%2F%2Fdoi.org%2F10.1038%252Fs41586-019-1272-6)](https://doi.org/10.1038%2Fs41586-019-1272-6)
- Cleavage Profile [![DOI](https://img.shields.io/badge/DOI-10.1073%2Fpnas.2209852119-blue?style=flat-square)](https://doi.org/10.1073/pnas.2209852119)

### Documentation
Documentation for FinaleToolkit can be found [here](https://epifluidlab.github.io/finaletoolkit-docs/).

### Wiki/Tutorials
The wiki and tutorial page for FinaleToolkit can be found [here](https://github.com/epifluidlab/FinaleToolkit/wiki).

### Compatible File Formats

FinaleToolkit is compatible with almost any paired-end sequence data:

- Binary Alignment Map (`.bam`) files with an associated index file (`.bam.bai`).
- Compressed Reference-oriented Alignment Map (`.cram`) files.
- Fragment (`.frag.gz`) files with an associated tabix index file (`.frag.gz.tbi`).

### Using Fragment Files

Fragment (`.frag.gz`) files are block-gzipped BED3+2 files with the following columns: `chrom` , `start` , `stop` , `mapq` , `strand`.

We encourage you to use our comprehensive database, FinaleDB, to access relevant fragment files. Learn more about FinaleDB [here](http://finaledb.research.cchmc.org).

## Contact
- James Li: lijw21@wfu.edu
- Ravi Bandaru: ravi.bandaru@northwestern.edu
- Yaping Liu: yaping@northwestern.edu

## License
This project falls under an MIT license. See the included `LICENSE` file for details.
