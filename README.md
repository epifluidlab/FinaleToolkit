# FinaleToolkit (rewrite)

**FinaleToolkit** (FragmentatIoN AnaLysis of cEll-free DNA Toolkit) is a package
and standalone program to extract fragmentation features of cell-free DNA from
paired-end sequencing data.

This is a full refactor of the original
[FinaleToolkit](https://github.com/epifluidlab/FinaleToolkit). The **Python API
and numeric behavior are preserved** (the original test suite passes unchanged),
while the code is modernized: a single import namespace, type hints, named
return types, an informative exception hierarchy, vectorized/streaming hot
paths, and a clear I/O / compute / CLI separation with shared helpers.

The **CLI subcommands are unchanged, but the flag names were redesigned** for
consistency (a clean break from the original — see `CHANGELOG.md` for the full
old→new flag mapping).

## Installation

```bash
pip install .
# or, for development:
pip install -e .
```

Requires Python ≥ 3.10. Key dependencies: numpy, pysam, pyBigWig, py2bit, numba,
scipy, pandas, statsmodels, loess, matplotlib, tqdm.

## Python API

Everything is reachable from the top-level namespace — no need to know internal
submodule paths:

```python
import finaletoolkit as ftk

# Coverage over intervals (returns NamedTuples that still unpack like the old tuples)
results = ftk.coverage("sample.bam", "intervals.bed", output_file=None)
for r in results:
    print(r.contig, r.start, r.stop, r.coverage)   # named access
    chrom, start, stop, name, cov = r               # tuple unpacking still works

# Windowed protection score
scores = ftk.wps("sample.bam", "chr12", 34443000, 34448000, chrom_size=133851895)

# End motifs + motif diversity score
motifs = ftk.end_motifs("sample.bam", "hg38.2bit", k=4)
print(motifs.motif_diversity_score())

# DELFI, cleavage profile, fragment lengths, breakpoint motifs, genome gaps ...
df = ftk.delfi("sample.bam", "autosomes.chrom.sizes", "bins.bed", "hg19.2bit", gap_file="hg19")
gaps = ftk.GenomeGaps.hg38()
```

The original import paths continue to work unchanged:

```python
from finaletoolkit.frag import wps, multi_wps, delfi, end_motifs, EndMotifFreqs
from finaletoolkit.utils import frag_generator, frag_array, agg_bw, filter_file
from finaletoolkit.genome import GenomeGaps, ContigGaps
from finaletoolkit.io import ReferenceWrapper, AlignmentWrapper, Fragment
```

## Command-line interface

```bash
finaletoolkit --help
finaletoolkit coverage sample.bam intervals.bed -o coverage.bed
finaletoolkit wps sample.bam tss.bed --chrom-sizes hg38.chrom.sizes -o wps.bw -t 8
finaletoolkit end-motifs sample.bam hg38.2bit -k 4 -o motifs.tsv
finaletoolkit delfi sample.bam autosomes.chrom.sizes hg19.2bit bins.bed -g hg19 -o delfi.tsv
```

Every subcommand has `--help` with an example invocation. The flag scheme is
consistent across all subcommands:

| Concept | Flag |
|---|---|
| Output path | `-o, --output` |
| Reference (optional, for CRAM) | `-r, --reference` |
| Minimum mapping quality | `-q, --min-mapq` |
| Fragment length bounds | `--min-length`, `--max-length` |
| Worker processes | `-t, --threads` |
| Verbosity (repeatable) | `-v, --verbose` |
| k-mer length | `-k, --kmer-length` |
| Strand selection (motifs) | `--strand {both,forward,reverse}` |
| Boolean on/off options | `--x` / `--no-x` (e.g. `--savgol`/`--no-savgol`) |

> **Note:** the command-line flag names are a clean break from the original
> FinaleToolkit. The **Python API is unchanged** — only CLI spellings changed.
> See `CHANGELOG.md` for the full old→new flag mapping.

### Subcommands

`coverage`, `frag-length-bins`, `frag-length-intervals`, `cleavage-profile`,
`wps`, `adjust-wps`, `delfi`, `delfi-gc-correct`, `end-motifs`,
`interval-end-motifs`, `breakpoint-motifs`, `interval-breakpoint-motifs`, `mds`,
`interval-mds`, `filter-file`, `agg-bw`, `gap-bed`.

## Supported input formats

- **Alignments:** indexed BAM and CRAM. CRAM requires a reference
  (`-r/--reference` on the CLI, or the `reference_file` argument in the API).
- **Fragments:** tabix-indexed FinaleDB `.frag.gz` (and BED6 `.bed.gz`, read
  with a warning).
- **Reference:** `.2bit` and FASTA (`.fa`/`.fasta`/`.fna`, optionally `.gz`;
  a `.fai` index is created automatically when missing).

## Package layout

```
finaletoolkit/
├── io/        # ReferenceWrapper, AlignmentWrapper/Fragment, output helpers
├── frag/      # feature extractors (frag length, coverage, WPS, cleavage, DELFI, motifs)
├── genome/    # GenomeGaps / ContigGaps + bundled UCSC gap tracks
├── utils/     # fragment streaming, intervals, k-mers, filtering, aggregation, parallelism
├── cli/       # argparse front-end (commands/ builders + dispatch)
└── exceptions.py  # FinaleToolkitError hierarchy
```

## License

MIT (see `LICENSE`).
