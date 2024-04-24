# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.0] - 2024-04-24

### Added
- Added finaletools.interval_end_motifs function to calculate end-motifs
over genomic intervals. Stores results in an IntervalEndMotifs object.
- Added CLI subcommand interval-end-motifs to calculate end-motifs over
genomic intervals.
- Added CLI subcommand interval-mds to calculate MDS over intervals from
interval end-motifs table.

### Changed
- Added gc_correct option to delfi_merge_bins so that merging is possible
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