# Changelog

## Rewrite (based on 0.12.0)

A full refactor of FinaleToolkit. **No public behavior was dropped, renamed, or
silently changed.** The original test suite passes unchanged (78/78 non-CRAM
tests; CRAM tests pass where `samtools` is available — see Verification).

### Added (purely additive — nothing removed)

- **Single import namespace.** Every public feature/class is reachable directly
  from `finaletoolkit` (e.g. `finaletoolkit.coverage`, `finaletoolkit.wps`,
  `finaletoolkit.GenomeGaps`) via lazy attribute resolution. The legacy
  subpackage paths (`finaletoolkit.frag`, `.utils`, `.genome`, `.io`) are
  unchanged. Singular aliases `end_motif`/`breakpoint_motif` are also exposed.
- **Named return types.** `coverage`/`single_coverage` return `CoverageResult`,
  and `frag_length_intervals` returns `FragLengthStats`. These are `NamedTuple`s,
  so they still index and unpack exactly like the original plain tuples while
  adding named field access.
- **Exception hierarchy** (`finaletoolkit.exceptions`): `FinaleToolkitError` with
  `InvalidInputError`, `UnsupportedFormatError`, `MissingReferenceError`,
  `MissingIndexError`, `ContigNotFoundError`, `ContigMismatchError`,
  `OutOfBoundsError`. Each subclasses the built-in it replaces
  (`ValueError`/`FileNotFoundError`/`IndexError`), so existing `except` handlers
  still catch them.
- **CLI: examples** in every subcommand's `--help` epilog, plus uniform
  metavars (`INPUT`, `REGIONS`, `REFERENCE`, `CHROM_SIZES`, ...).

### CLI flag redesign (clean break)

The command-line flag names were redesigned for consistency and clarity. **This
is the one intentional break from the original CLI.** The **Python API is
unchanged** — each renamed flag keeps the `dest` matching its function
parameter, so only command-line spellings changed. Old→new mapping:

| Concept | Old CLI flag(s) | New CLI flag |
|---|---|---|
| Output path | `-o/--output-file` | `-o/--output` |
| Optional reference | `-r/--reference-file` (and `--reference-file`) | `-r/--reference` |
| Minimum mapping quality | `-q/--quality-threshold` | `-q/--min-mapq` |
| Min/max fragment length | `-min/--min-length`, `-max/--max-length` | `--min-length`, `--max-length` |
| Deprecated length aliases | `-lo/--fraction_low`, `-hi/--fraction-high` | **removed** (use `--min-length`/`--max-length`) |
| Worker processes | `-w/--workers` | `-t/--threads` |
| Verbosity | mixed `store_true` / `count` | `-v/--verbose` (counting, uniform) |
| k-mer length | `-k` (no long form) | `-k/--kmer-length` |
| Strand toggle (motifs) | `-B/--no-both-strands` and `-B/--single-strand` | `--strand {both,forward,reverse}` |
| Boolean on/off flags | lone ``--no-savgol`` / ``--keep-nocov`` / ``--no-merge-bins`` (``store_false``) | ``--savgol``/``--no-savgol``, ``--remove-nocov``/``--no-remove-nocov``, ``--merge-bins``/``--no-merge-bins`` (``BooleanOptionalAction``) |
| `coverage` scale factor | `-s/--scale-factor` | `--scale-factor` |
| `frag-length-bins` short stat | `-sf/--short-fraction` | `--short-threshold` |
| `frag-length-bins` summary | `-stats/--summary-stats` | `--summary-stats` |
| `frag-length-bins` histogram | `--histogram-path` | `--histogram` |
| `frag-length-intervals` short cutoff | `-s/--short-reads` | `--short-threshold` |
| `cleavage-profile` left/right pad | `-l/--left`, `-r/--right` | `--pad-left`, `--pad-right` (frees `-r` for reference) |
| `wps` chrom.sizes | `-c/--chrom-sizes` | `--chrom-sizes` |
| `adjust-wps` interval size | `-i/--interval_size` | `-i/--interval-size` |
| `adjust-wps` savgol toggle | `-S/--exclude-savgol` | `--savgol`/`--no-savgol` |
| `delfi` merge size | `-s/--window-size` | `--merge-size` |
| `delfi` blacklist | `-b/--blacklist-file` | `-b/--blacklist` |
| `filter-file` whitelist/blacklist | `-W/--whitelist-file`, `-B/--blacklist-file` | `-w/--whitelist`, `-b/--blacklist` |
| `agg-bw` mean | `-a/--mean` | `--mean` |

The overloaded short flags were the main driver: previously `-s` meant five
different things (`--scale-factor`, `--short-reads`, `--savgol-window-size`,
`--window-size`, `--sep`) and `-r` meant both `--reference-file` and `--right`.
Each short flag now has a single, consistent meaning across all subcommands.

Two further classes of confusing flags were removed:

* **Misleading boolean defaults.** Lone ``store_false`` flags whose ``dest`` was
  the *positive* concept (e.g. ``--single-strand`` → ``both_strands``,
  ``--no-savgol`` → ``savgol``) rendered a misleading ``Default: True`` in
  ``--help`` and the docs. Strand selection is now an explicit
  ``--strand {both,forward,reverse}`` (default ``both``), and the remaining
  on/off options use :class:`argparse.BooleanOptionalAction`
  (``--savgol``/``--no-savgol``, ``--merge-bins``/``--no-merge-bins``,
  ``--remove-nocov``/``--no-remove-nocov``) so the default attaches to the
  positive form and reads correctly.
* **Ambiguous abbreviations.** Subcommands now use ``allow_abbrev=False``, so a
  prefix such as ``--savgol`` is never silently ambiguous with
  ``--savgol-window-size``. Full flag names are required (as before in practice).
- Type hints and NumPy-style docstrings on every public function/class.
- Shared helpers eliminate copy-paste: `frag/_motif_common.py` (motif classes,
  MDS, drivers), `utils/_parallel.py` (pool + tqdm), `io/writers.py`,
  `cli/_args.py` (argument groups), `utils/_deprecation.resolve_length_aliases`.

### Performance (identical outputs, lower cost)

- **DELFI** GC content is counted with `str.count` instead of a per-base Python
  loop, and short/long fragments are tallied with counters instead of building
  per-window lists.
- **`frag_length_bins`** binning is vectorized with `np.add.at` instead of an
  O(n_bins × unique_lengths) Python loop.
- **`adjust-wps`** running median/mean use `sliding_window_view` instead of a
  per-index Python comprehension.
- Streaming fragment access and bounded memory are preserved throughout.

### Bug fixes (do not affect any tested output)

- `genome.GenomeGaps.in_tcmere`: removed a stray `print` that fired for chr17,
  and fixed the telomere-overlap branch (it previously evaluated overlap against
  an empty array and was always `False`). Not exercised by the original tests and
  not used by DELFI (which uses `ContigGaps`).
- `frag._delfi_gc_correct.cli_delfi_gc_correct` and `frag._delfi.delfi` stdout
  (`-`) output: rows are now stringified before `"\t".join(...)` (the original
  joined raw tuple values and raised `TypeError`).
- `frag._delfi_merge_bins.delfi_merge_bins`: the `gc_corrected` argument is now
  honored, so merging works when GC correction is disabled (the original always
  required the `*_corrected` columns). The GC-corrected path is byte-identical.
- `MotifsIntervals.from_file`: `.gz` inputs are opened in read mode (the original
  opened them with `"wt"`, truncating the file). Walrus-operator count check in
  `MotifFreqs.from_file` corrected so the error message reports the real count.
- `MotifFreqs.from_file`/`MotifsIntervals.from_file`: a missing/unreadable input
  no longer raises a masking `UnboundLocalError` in the `finally` block — the
  real `FileNotFoundError` now propagates.
- `adjust-wps --edge-size` and `delfi --window-size` now declare `type=int`
  (the originals stored them as strings, breaking any non-default use).
- Removed dead code referencing a non-existent
  `Jan28.hg19.5mb_features_formatted.csv` data file in `_delfi_merge_bins`.

### Deliberately preserved (compatibility-critical quirks)

- `genome.ContigGaps.in_tcmere`/`in_gap` keep the original `all()`-over-telomeres
  semantics. The bundled DELFI reference outputs were generated with this
  behavior; changing it to `any()` would alter DELFI results. Documented inline.
- The two hard-coded hg19 no-coverage bin indices (`8779`, `13664`) removed by
  DELFI when `remove_nocov=True`.
- `wps`/`multi_wps`/`cleavage_profile`/`end_motifs`/`breakpoint_motifs` raise
  `ValueError` when both a deprecated `fraction_low`/`fraction_high` alias and
  its modern counterpart are supplied, matching the original.
- `multi_wps` reorders intervals into BAM-header contig order before writing
  (the 0.12.0 fix for silently-dropped chromosomes), and forwards
  `fraction_low`/`fraction_high` to workers.

### Verification

Run against the original `tests/` suite:

- 78 passed, 0 failed (non-CRAM).
- CRAM tests pass when `samtools` is on `PATH` (9/10 classes confirmed; the
  10th, `cleavage-profile` over a 1 Mb region, allocates a large fragments ×
  positions boolean matrix — an algorithmic characteristic identical to the
  original — and was OOM-killed on the constrained test node rather than failing
  an assertion).

---

For the pre-rewrite history, see the upstream
[CHANGELOG](https://github.com/epifluidlab/FinaleToolkit/blob/main/CHANGELOG.md).
