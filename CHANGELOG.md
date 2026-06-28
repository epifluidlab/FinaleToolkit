# Changelog

All notable changes to this project will be documented in this file.

The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-06-26

A full refactor and modernization of FinaleToolkit. The **Python API and numeric
behavior are preserved** — the test suite passes (80 passed, 27 skipped; CRAM
tests run where `samtools` is available). The internals, command line, and
documentation were modernized.

### Added

- **Single import namespace.** Every public feature/class is reachable directly
  from `finaletoolkit` (e.g. `finaletoolkit.coverage`, `finaletoolkit.wps`,
  `finaletoolkit.GenomeGaps`) via lazy attribute resolution. The subpackage
  paths (`finaletoolkit.frag`, `.utils`, `.genome`, `.io`) still work. Singular
  aliases `end_motif`/`breakpoint_motif` are also exposed.
- **Named return types.** `coverage`/`single_coverage` return `CoverageResult`
  and `frag_length_intervals` returns `FragLengthStats` — `NamedTuple`s that
  index and unpack like plain tuples while adding documented named fields.
- **Exception hierarchy** (`finaletoolkit.exceptions`): `FinaleToolkitError`
  with `InvalidInputError`, `UnsupportedFormatError`, `MissingReferenceError`,
  `MissingIndexError`, `ContigNotFoundError`, `ContigMismatchError`,
  `OutOfBoundsError`. Each subclasses the built-in it replaces
  (`ValueError`/`FileNotFoundError`/`IndexError`), so existing `except` handlers
  still catch them.
- **Click command line, rendered with rich-click.** Every subcommand's `--help`
  shows panels grouping related options, an example invocation, and uniform
  metavars (`INPUT`, `REGIONS`, `REFERENCE`, `CHROM_SIZES`, ...). Accent colors
  are chosen to stay legible on both dark and light terminals. Pass `-o -` to
  write any output to standard output (stdout).
- **Type hints and NumPy-style docstrings** on every public function/class.
- Shared helpers eliminate copy-paste: `frag/_motif_common.py` (motif classes,
  MDS, drivers), `utils/_parallel.py` (pool + tqdm), `io/writers.py`,
  `cli/_args.py` (reusable Click options), `cli/_dispatch.py` (lazy dispatch).

### Changed

- **Command line reimplemented on Click + rich-click** (previously argparse).
  The flag names were redesigned for consistency — the one intentional break
  from the previous CLI. The **Python API is unchanged**: each flag keeps the
  parameter name matching its function argument, so only command-line spellings
  changed.

  | Concept | Old CLI flag(s) | New CLI flag |
  |---|---|---|
  | Output path | `-o/--output-file` | `-o/--output` (`-` = stdout) |
  | Optional reference | `-r/--reference-file` | `-r/--reference` |
  | Minimum mapping quality | `-q/--quality-threshold` | `-q/--min-mapq` |
  | Min/max fragment length | `-min/--min-length`, `-max/--max-length` | `--min-length`, `--max-length` |
  | Deprecated length aliases | `-lo/--fraction_low`, `-hi/--fraction-high` | **removed** |
  | Worker processes | `-w/--workers` | `-t/--threads` |
  | Verbosity | mixed `store_true`/`count` | `-v/--verbose` (counting) |
  | k-mer length | `-k` (no long form) | `-k/--kmer-length` |
  | Strand toggle (motifs) | `-B/--no-both-strands`, `-B/--single-strand` | `--strand {both,forward,reverse}` |
  | `coverage` scale factor | `-s/--scale-factor` | `--scale-factor` |
  | `frag-length-bins` short stat | `-sf/--short-fraction` | `--short-threshold` |
  | `frag-length-bins` summary | `-stats/--summary-stats` | `--summary-stats` |
  | `frag-length-bins` histogram | `--histogram-path` | `--histogram` |
  | `frag-length-intervals` short cutoff | `-s/--short-reads` | `--short-threshold` |
  | `cleavage-profile` left/right pad | `-l/--left`, `-r/--right` | `--pad-left`, `--pad-right` |
  | `wps` chrom.sizes | `-c/--chrom-sizes` | `--chrom-sizes` |
  | `adjust-wps` interval size | `-i/--interval_size` | `-i/--interval-size` |
  | `adjust-wps` savgol toggle | `-S/--exclude-savgol` | `--savgol/--no-savgol` |
  | `delfi` merge size | `-s/--window-size` | `--merge-size` |
  | `delfi` blacklist | `-b/--blacklist-file` | `-b/--blacklist` |
  | `filter-file` whitelist/blacklist | `-W/--whitelist-file`, `-B/--blacklist-file` | `-w/--whitelist`, `-b/--blacklist` |
  | `agg-bw` mean | `-a/--mean` | `--mean` |

  The overloaded short flags were the driver: `-s` previously meant five
  different things and `-r` meant two; each short flag now has one meaning.
  Boolean options use Click's `--x/--no-x` pairs (e.g. `--savgol/--no-savgol`)
  with the default stated in the help text, and subcommands disable prefix
  matching so `--savgol` is never ambiguous with `--savgol-window-size`.
- **`interval-mds` renamed to `regional-mds`** — the regional Motif Diversity
  Score (rMDS), per Bandaru et al., *Journal of Clinical Investigation* 2026
  ([196284](https://www.jci.org/articles/view/196284)). The implementing
  function is `_cli_regional_mds`.
- **Documentation rebuilt** on a custom, compact in-tree Sphinx theme
  (Northwestern-purple accents, Inter typography, light/dark). The CLI reference
  is generated from the live Click commands (sphinx-click), the API reference
  from docstrings (autodoc), and the user guide was restructured. CI
  (`build-docs.yml`, `.readthedocs.yaml`) installs the package and the docs
  toolchain.

### Removed

- The deprecated `delfi-gc-correct` CLI command. GC correction is performed by
  `delfi` automatically (`delfi --no-gc-correct` opts out). The Python function
  `finaletoolkit.frag.delfi_gc_correct` is unchanged and still available.

### Performance (identical outputs, lower cost)

- **DELFI** worker pool is ~80x faster on whole-genome inputs while producing
  bit-identical output: the blacklist BED is parsed once and filtered per window
  with binary search (instead of being re-read and linearly scanned for every
  100kb window), and one alignment handle and one reference handle are opened
  per worker via the `Pool` initializer and reused for every window. Per-contig
  `ContigGaps` are preloaded into worker globals rather than pickled into every
  task. Contributed by D.H.K. (Duco) Gaillard
  ([@DucoG](https://github.com/DucoG)) in
  [#172](https://github.com/epifluidlab/FinaleToolkit/pull/172); guarded by new
  `test_workers_equivalence` and `test_fragfile_input` tests.
- **DELFI** GC content is counted with `str.count` instead of a per-base Python
  loop, and short/long fragments are tallied with counters.
- **`frag_length_bins`** binning is vectorized with `np.add.at`.
- **`adjust-wps`** running median/mean use `sliding_window_view`.
- Streaming fragment access and bounded memory are preserved throughout.

### Fixed

- `genome.GenomeGaps.in_tcmere`: removed a stray `print` that fired for chr17,
  and fixed the telomere-overlap branch (it previously evaluated overlap against
  an empty array and was always `False`).
- `frag._delfi_gc_correct.cli_delfi_gc_correct` and `frag._delfi.delfi` stdout
  (`-`) output: rows are stringified before `"\t".join(...)` (previously joined
  raw tuple values and raised `TypeError`).
- `frag._delfi_merge_bins.delfi_merge_bins`: the `gc_corrected` argument is
  honored, so merging works when GC correction is disabled.
- `MotifsIntervals.from_file`: `.gz` inputs are opened in read mode (previously
  opened with `"wt"`, truncating the file). The walrus-operator count check in
  `MotifFreqs.from_file` reports the real count.
- `MotifFreqs.from_file`/`MotifsIntervals.from_file`: a missing/unreadable input
  no longer raises a masking `UnboundLocalError`; the real `FileNotFoundError`
  propagates.
- `adjust-wps --edge-size` and `delfi --window-size` declare `type=int`.

### Preserved (compatibility-critical quirks)

- `genome.ContigGaps.in_tcmere`/`in_gap` keep the `all()`-over-telomeres
  semantics that the bundled DELFI reference outputs were generated with.
- The two hard-coded hg19 no-coverage bin indices (`8779`, `13664`) removed by
  DELFI when `remove_nocov=True`.
- `multi_wps` reorders intervals into BAM-header contig order before writing
  (the 0.12.0 fix for silently-dropped chromosomes) and forwards
  `fraction_low`/`fraction_high` to workers.

## [0.12.0] - 2026-05-28

### Added
- CRAM support threaded through all `frag` modules: `coverage`, `frag_length`,
  `frag_length_bins`, `frag_length_intervals`, `wps`, `multi_wps`,
  `cleavage_profile`, and `multi_cleavage_profile` now accept a
  `reference_file` parameter (FASTA only) required for decoding CRAM input.
- `--reference-file` CLI flag added to `frag-length-bins`,
  `frag-length-intervals`, and `wps` subcommands.
- `--reference-file` CLI flag added to `cleavage-profile` subcommand (no `-r`
  short form, as `-r` was already taken by `--right`).
- `ReferenceWrapper` auto-creates a missing `.fai` index (via `pysam.faidx`)
  when a FASTA file is opened without a pre-built index.
- `.fna` and `.fna.gz` extensions recognised as FASTA by `ReferenceWrapper`.
- Integration test suite `tests/test_cram.py` verifying CRAM produces
  identical results to BAM for `single_coverage`, `frag_length_bins`, `wps`,
  and `delfi`.
- Smoke tests in `tests/test_cram.py` verifying CRAM input runs without error
  for `coverage`, `frag_length`, `frag_length_intervals`, `multi_wps`,
  `cleavage_profile`, and `multi_cleavage_profile`.

### Changed
- CLI help text for `delfi`, `end-motifs`, `interval-end-motifs`,
  `breakpoint-motifs`, and `interval-breakpoint-motifs` updated: the
  `reference_file` argument description now reads "A .2bit or FASTA
  (.fa, .fasta, .fna) file" instead of "The .2bit file".
- `multi_wps` now correctly forwards `fraction_low` and `fraction_high` to
  worker processes (these were silently dropped before).

### Fixed
- CRAM files can now be used with all analysis subcommands; previously
  `frag_generator` was never passed the reference file needed for CRAM
  decoding.
- `delfi` docstring clarifies that when `input_file` is a CRAM file,
  `reference_file` must be a FASTA (not .2bit), as htslib requires FASTA for
  CRAM decoding.
- `multi_wps` (and therefore `wps` CLI) no longer silently drops chromosomes
  from BigWig output when the input BED file is sorted in a different
  chromosome order than the BAM header (e.g. alphabetical vs. numeric). The
  interval list is now sorted by BAM-header chromosome order before writing,
  preventing the `RuntimeError` that previously caused silent data loss for
  chr2–chr9, chrX, chrY and similar chromosomes.

## [0.11.1] - 2026-04-21

### Added
- Python 3.13 added to GitHub Actions CI test matrix.
- Snakemake workflow section added to README.

### Fixed
- `multi_wps` no longer crashes with `ValueError: negative dimensions are not
  allowed` when BED intervals produce degenerate windows (e.g. when expanded
  intervals exceed chromosome boundaries or when closely spaced TSS sites
  trigger overlap-trimming that inverts `start`/`stop`).
  - `stop` is now clamped to chromosome size when computing expanded windows.
  - Degenerate intervals (`stop <= start`) produced by overlap-trimming are
    silently skipped instead of being forwarded to worker processes.
  - `wps` now checks for `stop <= start` at function entry, emits a
    `UserWarning`, and returns an empty result rather than crashing.

### Changed
- `low_quality_read_pairs` in `utils` now also checks the `MQ` (mate mapping
  quality) BAM tag and marks a read pair as low quality if the mate's mapping
  quality falls below `min_mapq`.

## [0.11.0] - 2025-09-05

### Added
- New `breakpoint-motifs` and `interval-breakpoint-motifs` CLI subcommands and
  corresponding Python API (`finaletoolkit.frag.breakpoint_motifs`) for
  computing breakpoint motif frequencies from cfDNA fragments.
- `--summary-stats` flag for `frag-length-bins`: appends summary statistics
  (mean, median, short fraction, etc.) as comment lines to the output file.
- `--short-fraction` option for `frag-length-bins` to specify the length
  threshold used when computing the short-fragment fraction.
- `--short-reads` option for `frag-length-intervals` to customize the short
  read length threshold (previously hard-coded).
- Automated documentation build workflow (`.github/workflows/build-docs.yml`).
- `frag.gz` (FinaleDB format) accepted wherever `bed.gz` fragment files were
  previously accepted in `frag_generator`.

### Changed
- `chrom_sizes` is now a **required** argument for Python API functions and CLI
  commands that previously accepted it as optional (affects `multi_wps`,
  `cleavage_profile`, and related subcommands).
- DELFI GC-correction flag logic made more intuitive: passing `-G` /
  `--no-gc-correct` now unambiguously disables GC correction.
- Deprecated DELFI subcommands now appear in verbose logs.

### Fixed
- `cleavage-profile` now correctly reads chromosome sizes when a `chrom_sizes`
  file is provided (PR #158).
- `multi_wps` now issues a `UserWarning` and skips BED intervals whose
  chromosome is absent from `chrom_sizes`, rather than raising a `KeyError`.
- Additional error checking for intervals passed to `wps` (e.g. `start >
  stop` now raises a descriptive error).
- `end_motifs` no longer duplicates the last interval when parallelizing over
  chromosomes.
- `breakpoint_motifs` parallelization race condition fixed.
- `delfi` and `adjust-wps` now handle chromosomes without centromere
  annotations (e.g. chrM, chrX, chrY) without crashing (PR #156, @skchronicles).
- `agg_bigwig` no longer crashes when a BigWig file has no entries over a
  given genomic range (PR #156, @skchronicles).

## [0.10.7] - 2025-1-22

### Added
- Snakemake workflow included in FinaleToolkit documentation.

### Fixed
- Issue with `filter-file` where the lack of `output_file.flush()` would
create incorrect outputs.

## [0.10.6] - 2025-1-14

### Fixed
- More descriptive error message in `multi_wps` when an invalid interval
with `start > stop` is encountered. Now the chromosome name is mentioned
to assist users.

## [0.10.5] - 2025-1-9

### Changed
- added `intersect_policy` to `filter_file`

### Fixed
- blacklist functionality of `filter_file`

## [0.10.4] - 2025-1-2

### Changed
- Changed `filter-bam` function into `filter-file`; Now accepts BED and CRAM in
addition to BAM
- filter-file command can now accept a blacklist file

## [0.10.3] - 2024-12-21

### Fixed
- Update various docstrings and help statements to no longer mention SAM
files as an accepted format

## [0.10.2] - 2024-12-19

### Fixed
- changed default args for `wps` would lead to errors. Now `wps` defaults to
LWPS fragment lengths (120-180nt).

### Changed
- made `finaeltoolkit.utils.typing` public. This is a module containing some
useful type aliases
- minor formatting and typing changes
- renamed `fraction_low` and `fraction_high` to `min_length` and `max_length`
for `.utils.frag_array` and `.frag.wps`. `wps` retains the deprecated arg names
but issues a warning.

## [0.10.1] - 2024-12-19

### Fixed
- Added missing `-n` arg to `end-motifs`.
- Fixed incorrect `ValueError` regarding the `negative_strand` arg.
- Incorrect function name for `wps` leading to errors when called from CLI.

### Added
- Additional tests for the CLI lazy loading implementation

## [0.10.0] - 2024-12-18

### Changed
- several modules containing implementations of fragmentomic features or
utiliy functions have been made internal. This means there is now only one
obvious import for each function. For example, `multi_wps` is imported from
`finaletoolkit.frag`, and no longer can be imported from
`finaletoolkit.frag.multi_wps`
- The CLI now uses lazy importing, drastically speeding up finaletoolkit when
called from a command line.
- Added `negative_strand` option for end motifs related functions. When
used in conjunction with `both_strands`, only end motifs on the negative
(Crick) strand are considered in calculations.
- Renamed `fraction_high` and `fraction_low` in `utils.utils.frag_generator`
to `min_length` and `max_length`.

### Fixed
- deprecated arguments for `end-motifs` had default values which could
lead to an error. This is fixed.

## [0.9.1] - 2024-12-17

### Fixed
- CLI no longer prints an error message if `finaletoolkit` is called without args.
- `frag-length-bins`, when writing a file, now writes the interval between
`min` and `max` as inclusive. That is, previously when `min=1` and `max=2`,
only fragments of length 1 are reported. Now when such a result is calculated, 
the interval given is `min=1` and `max=1`.
- Updated some descriptions and docstrings.

### Added
- `adjust-wps` now has an option `-S` or `--exclude-savgol` to not perform
Savitsky-Golay filtering.

### Changed
- Several CLI options were renamed so that underscores become hyphens. This is
for consistency and to simplify writing commands.

## [0.9.0] - 2024-12-16

### Removed
- `strand_location` arg from `agg_bigwig`
- `cli_hist` module

### Fixed
- fixed bug involving tqdm progress bar in `frag_length_intervals`
- some code formatting
- fixed bug involving arg names in `filter-bam`
- add some missing args to CLI
- issues with running `cleavage-profile` (#115)
- issues with writing to bigwig with `wps`

### Changed
- change default of arg `both_strands` of `end_motifs` to True to match
behavior of original scripts
- rename `fraction_high` and `fraction_low` to  `min_length` and `max_length`
for all features, deprecating old args as aliases if needed.
- numpy 2 compatible
- fragmentomics functions assume Tabix indexed files all follow the
FinaleDB Frag.gz file format. That is, columns are `chrom`, `start`, `stop`,
`score`, and `strand`. If more columns are detect, a warning is issued, and
FinaleToolkit will attempt to parse the file as a BED6 format.
- renamed `genome_file` to `chrom_sizes` for most functions.
- `multi_wps` and `multi_cleavage_profile` no longer return a value due to memory issues when attempting to calculate these genomewide. Instead, users should refer to the file specified with `output_file`.

### Added
- internal `utils._typing` and `utils._deprecation` modules
- test for `delfi`

### Deprecated
- `delfi-gc-correct` command. GC-correction is performed automatically by `delfi`
already.

## [0.8.0] - 2024-12-04

### Removed
- `finaletoolkit.frag.frag_length_bins` no longer has the `contig_by_contig`
option. This never had any functionality.
- `finaletoolkit.frag.frag_length_bins` no longer generates a text-based
histogram.

### Fixed
- `contig_sizes` option included for `cleavage-profile` CLI command.
- `normalize` option for `coverage` fixed so it no longer normalizes twice
- `normalize=False` for `coverage` runs much faster
- misc typehints and docstrings

### Changed
- `finaletoolkit.frag.frag_length_bins` uses a dict based implementation
that is more memory efficient.
- `finaletoolkit.frag.frag_length_bins` and
`finaletoolkit.frag.frag_length_intervals` now take `min_length` and 
`max_length` keyword args to only consider fragments of certain lengths.
- flags for `frag-length-bins` and `frag-length-intervals` CLI commands updated
to match Python API
- `coverage` default argument for `normalize` changed to `False` 
- `coverage` default argument for `scale_factor` changed to 1. 

### Added
- `finaletoolkit.frag.frag_length_bins` can generate a histogram figure
- tests for `frag_length` module

## [0.7.8] - 2024-11-28

### Fixed
- update docs, docstring, and help message for wps to mention that
`site_bed` must be sorted.

### Added
- `normalize` keyword argument and `--normalize` flag to
`finaletoolkit.frag.coverage` function and `finaletoolkit coverage` subcommand,
respectively. Setting this argument/flag to true results in the output
being normalized by the total coverage, ignoring `scale_factor` if specified.
- `--intersect-policy` or `-p` flag added to `finaletoolkit coverage`subcommand.

## [0.7.7] - 2024-11-27

### Fixed
- subpackages can now be accessed when importing `finaletoolkit`. Previously,
the following code resulted in an error:
```python
>>> import finaletoolkit as ftk
>>> help(ftk.frag)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
AttributeError: module 'finaletoolkit' has no attribute 'frag'
```
Now this is a valid way to access subpackages `cli`, `frag`, `genome`, and
`utils`.

## [0.7.6] - 2024-11-18

### Fixed
- indexing issue in region_end_motifs that would misread strand
information when calculating end motifs on forward-strand only.

### Changed
- frag_generator now accepts fragment coordinates in bed.gz files

## [0.7.5] - 2024-10-10

### Fixed
 - `delfi` accepts `gap_file=None`
 - update prog for `delfi` to reflect compatibility with reference genomes
 other than hg19

## [0.7.4] - 2024-08-24

### Changed
 - Added many tests for util functions

### Fixed
 - Changed a nopython function to use numba compatible indexing

## [0.7.3] - 2024-08-20

### Changed
 - Used "not" instead of "~" in an if statement
 - Added a test for the coverage function

### Fixed
 - Ensured that the coverage value returns the expected value (previously
 returned an empty generator)

## [0.7.2] - 2024-08-17

### Changed
 - Included `output_file` as required argument for
 `finaletoolkit cleavage-profile`.

### Fixed
 - Fixed incompatible types in min function through an explicit cast of
 chrom_sizes to integers.

## [0.7.1] - 2024-08-11

### Changed
- include `chrom_sizes` file as required argument for
`finaletoolkit cleavage-profile`
- Numpy dependency version set to <2 to avoid breaking changes from numpy 2.
This will change in the future as we migrate to use numpy 2.
- Replaced all instances of `np.NaN` with `np.nan`.

### Fixed
- Fixed minor issues with typing in `finaletoolkit.genome.gaps`
- Fixed issue where data files are not packaged with FinaleToolkit

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
- `interval-mds` CLI subcommand calculates correctly without large negative
values.
- `interval-mds` CLI subcommand now correctly parses tsv files.

### Added
- Most end-motif related Python functions accept Path instances as inputs for
files.
- Unit and function tests, especially for end-motif related functions.

## [0.5.1] - 2024-05-03

### Changed
- All instances of finaletools have been renamed to finaletoolkit
- All default tabular files are now TSV
- Update contacts in TOML

### Fixed
- `interval-mds` and `mds` both calculate correctly when one motif has a
frequency of 0

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
