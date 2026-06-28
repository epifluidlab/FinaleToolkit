# <img alt="FinaleToolkit logo" src="docs/_static/finaletoolkit_logo_rounded.png" height="60"> &nbsp;FinaleToolkit

[![Docs](https://img.shields.io/badge/docs-epifluidlab.github.io-836eaa?style=flat-square)](https://epifluidlab.github.io/FinaleToolkit/)
[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fbioadv%2Fvbaf236-836eaa?style=flat-square)](https://doi.org/10.1093/bioadv/vbaf236)
[![License: MIT](https://img.shields.io/badge/License-MIT-836eaa?style=flat-square)](LICENSE)

<details>
<summary><b>Table of Contents</b></summary>

- [About The Project](#about-the-project)
- [Installation](#installation)
- [Usage](#usage)
  - [Functionality](#functionality)
  - [Documentation](#documentation)
  - [Wiki / Tutorials](#wiki--tutorials)
  - [Compatible File Formats](#compatible-file-formats)
  - [Using Fragment Files](#using-fragment-files)
  - [Snakemake Workflow](#snakemake-workflow)
- [Contact](#contact)
- [License](#license)

</details>

## About The Project

FinaleToolkit (**F**ragmentat**I**o**N** **A**na**L**ysis of c**E**ll-free DNA
**Toolkit**) is a package and standalone program to extract fragmentation
features of cell-free DNA from paired-end sequencing data.

### Citation

If you use FinaleToolkit in your research, please consider citing our paper:

> James Wenhan Li, Ravi Bandaru, Kundan Baliga, Yaping Liu.
> **FinaleToolkit: Accelerating Cell-Free DNA Fragmentation Analysis with a
> High-Speed Computational Toolkit.** *Bioinformatics Advances*, 2025, vbaf236.
> <https://doi.org/10.1093/bioadv/vbaf236>

## Installation

Install with `conda`:

```bash
conda install -c bioconda -c conda-forge finaletoolkit
```

Or with `pip`:

```bash
pip install finaletoolkit
```

## Usage

### Functionality

FinaleToolkit supports the following cell-free DNA fragmentation features:

- Fragment Length
- Coverage
- End Motifs [![DOI](https://img.shields.io/badge/DOI-10.1158%2F2159--8290.CD--19--0622-836eaa?style=flat-square)](https://doi.org/10.1158/2159-8290.cd-19-0622)
- Motif Diversity Score [![DOI](https://img.shields.io/badge/DOI-10.1158%2F2159--8290.CD--19--0622-836eaa?style=flat-square)](https://doi.org/10.1158/2159-8290.CD-19-0622)
- Windowed Protection Score [![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.cell.2015.11.050-836eaa?style=flat-square)](https://doi.org/10.1016/j.cell.2015.11.050)
- DELFI [![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41586--019--1272--6-836eaa?style=flat-square)](https://doi.org/10.1038/s41586-019-1272-6)
- Cleavage Profile [![DOI](https://img.shields.io/badge/DOI-10.1073%2Fpnas.2209852119-836eaa?style=flat-square)](https://doi.org/10.1073/pnas.2209852119)

### Documentation

Documentation for FinaleToolkit can be found
[here](https://epifluidlab.github.io/FinaleToolkit/).

### Wiki / Tutorials

The wiki and tutorial page for FinaleToolkit can be found
[here](https://github.com/epifluidlab/FinaleToolkit/wiki).

### Compatible File Formats

FinaleToolkit is compatible with almost any paired-end sequence data:

- Binary Alignment Map (`.bam`) files with an associated index file (`.bam.bai`).
- Compressed Reference-oriented Alignment Map (`.cram`) files.
- Fragment (`.frag.gz`) files with an associated tabix index file (`.frag.gz.tbi`).

### Using Fragment Files

Fragment (`.frag.gz`) files are block-gzipped BED3+2 files with the following
columns: `chrom`, `start`, `stop`, `mapq`, `strand`.

We encourage you to use our comprehensive database, FinaleDB, to access relevant
fragment files. Learn more about FinaleDB [here](http://finaledb.research.cchmc.org).

### Snakemake Workflow

Check out our
[Snakemake workflow](https://github.com/epifluidlab/finaletoolkit_workflow)!

## Contact

- James Li: james.li3@northwestern.edu
- Ravi Bandaru: ravi.bandaru@northwestern.edu
- Yaping Liu: yaping@northwestern.edu

## License

This project falls under an MIT license. See the included `LICENSE` file for
details.
