# Changelog

All notable changes to this project will be documented in this file.

The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased version]

### Changed
- include `chrom_sizes` file as required argument for
`finaletoolkit cleavage-profile`
- Numpy dependency version set to <2 to avoid breaking changes from numpy 2.0
- Replaced all instances of `np.NaN` with `np.nan`.

## [0.7.0] - 2024-07-21

### Added
- Brief description of modules in documentation under structure page
- Docstring
- `finaletoolkit.version` module containing single-source `__version__`
variable
- `remove_nocov` option in `finaletoolkit.frag.delfi` to toggle dropping two
bins with low coverage. These bins are dropped in delfi_scripts but
may not apply to fragment files not aligned to hg19.
- tests for `finaletoolkit.frag.delfi_merge_bins`

### Changed
- `finaletoolkit.frag.delfi` changed to accept files aligned to almost
any reference genome.
- `finaletoolkit.frag.delfi_merge_bins` algorithm changed to be
reference genome-agnostic and consistent with delfi_scripts
- `finaletoolkit delfi` options `-G`, `-M`, and `-R` to drop
gc-correction, merging, and remove nocov bins, respectively.

### Removed
- unused flags for `finaletools delfi`: `-W`, `--window-size`
- redundant flags for `finaletools delfi`: `-gc`, `--gc-correct`,
`-m`, `--merge-bins`

## [0.6.5] - 2024-07-15

### Changed
- `utils.agg_bw` now supports `PathLike` for input
- docstrings for `frag.end_motifs.EndMotifsIntervals` changed to be compatible
with Sphinx 

### Fixed
- added missing `gzip` import for `utils.agg_bw` 

### Added
- tests for `utils.agg_bw` 

## [0.6.4] - 2024-06-04

### Added
- `interval_size` argument for `adjust_wps`

### Changed
- `adjust_wps` checks if `median_window` is larger than interval
- remove default options from some private helper functions for better
error catching/predictable behavior.

### Fixed
- `wps`related functions and subcommands

## [0.6.3] - 2024-05-31

### Fixed
- fixed writing to `bed.gz` files when using `coverage`

## [0.6.2] - 2024-05-30

### Fixed
- adjusted handling of contig, start, stop for `frag_generator` so that
`coverage` does not throw exceptions for genomewide intervals.

### Added
- test for `single_coverage`

## [0.6.1] - 2024-05-26

### Changed
- add `__version` attribute
- `finaletoolkit --version` displays package version
- update PyPI page to include links

## [0.6.0] - 2024-05-26

### Fixed
- Fixed intersect policy for `cleavage_profile`. Now it calls `frag_generator`
with a policy of `any`.
- Clean up some comments and docstrings
- Fixed logging from coverage function

### Added
- Added numerous util functions
- Added `left` and `right` options to `cleavage_profile` and CLI
`cleavage-profile`.
- Added tests for cleavage profile and WPS.

### Changed
- Minimum Python version 3.9
- Changed `filter_bam` to have same filters as FinaleDB
- `utils.frag_generator` raises `ValueError` if `start` or `stop`
are specified without `contig`
- Type hints changed to use literals when possible

### Removed
- Removed `utils.get_contig_lengths`
- Removed `data`, `conda_envs`, and `figs` directories
- Removed unused dependencies `click`, `pybedtools`, and `cython
- Remove some unused imports from module files


## [0.5.2] - 2024-05-08

### Fixed
- `interval-mds` CLI subcommand calculates correctly without large negative values.
- `interval-mds` CLI subcommand now correctly parses tsv files.

### Added
- Most end-motif related Python functions accept Path instances as inputs for files.
- Unit and function tests, especially for end-motif related functions.


## [0.5.1] - 2024-05-03

### Changed
- All instances of finaletools have been renamed to finaletoolkit
- All default tabular files are now TSV
- Update contacts in TOML

### Fixed
- `interval-mds` and `mds` both calculate correctly when one motif has a frequency of 0

## [0.5.0] - 2024-04-24

### Added
- Added `finaletools.interval_end_motifs` function to calculate end-motifs
over genomic intervals. Stores results in an IntervalEndMotifs object.
- Added CLI subcommand `interval-end-motifs` to calculate end-motifs over
genomic intervals.
- Added CLI subcommand `interval-mds` to calculate MDS over intervals from
interval end-motifs table.

### Changed
- Added `gc_correct` option to `delfi_merge_bins` so that merging is possible
without GC correction

### Fixed
- `delfi` can now be run with `gc_correct=false` and `merge_bins=true`
- fixed `cleavage_profile` import in `frag`

## [0.4.5] - 2024-04-9

### Added
- Added `CHANGELOG.md`

### Fixed
- Fixed bug in coverage where  writing to non-bedgraph files would result in an
error

## [0.4.4] - 2024-04-5

### Changed
- `finaletools.frag.coverage` accepts Frag.gz format files
- update CLI help messages and docstrings for coverage and DELFI to reflect
current and previous changes
- update docs

## [0.4.2] - 2024-03-28

### Changed
- Updated emails in `pyproject.toml`
