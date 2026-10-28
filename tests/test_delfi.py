"""
Tests for modules related to calculating DELFI, including 
finaletoolkit.frag.delfi_merge_bins, finaletoolkit.frag.delfi, and
finaletoolkit.frag.delfi_gc_correct.
"""

import numpy as np
import pytest

import pandas as pd
import pysam
from click.testing import CliRunner

from finaletoolkit.cli.main_cli import main_cli
from finaletoolkit.frag import delfi_merge_bins
from finaletoolkit.frag import delfi
from finaletoolkit.frag._delfi import _resolve_gaps
from finaletoolkit.frag._delfi_merge_bins import _resolve_bin_size
from finaletoolkit.genome.gaps import GenomeGaps
from finaletoolkit.utils import frag_generator


def test_merge_bins(request):
    # getting files for comparison
    delfi_bins_csv = request.path.parent / 'data' / 'delfi' / 'test_delfi_100kb.csv'
    delfi_merged_bins_csv = request.path.parent / 'data' / 'delfi' / 'test_delfi_5mb.csv'

    delfi_bins = pd.read_csv(
        delfi_bins_csv, dtype={'contig':str, 'start':int, 'stop':int})
    delfi_merged_bins = pd.read_csv(
        delfi_merged_bins_csv, dtype={'contig':str, 'start':int, 'stop':int})

    merged_bins = delfi_merge_bins(delfi_bins)

    # same number of bins
    assert merged_bins.shape == delfi_merged_bins.shape

    # same bins
    assert (merged_bins['start'] == delfi_merged_bins['start']).all()
    assert (merged_bins['stop'] == delfi_merged_bins['stop']).all()

    # similar ratios
    assert (pytest.approx(merged_bins['ratio_corrected'], rel=5e-2)
            == delfi_merged_bins['ratio_corrected'])

def test_overall(request):
    """Testing entire delfi method on first 5Mb bin of hg19."""
    # getting files
    frag_file = request.path.parent / 'data' / 'delfi' / 'hg19.chr1.6Mb.bam'
    autosomes = request.path.parent / 'data' / 'delfi' / 'human.hg19.chr1.6Mb.genome'
    bins_file = request.path.parent / 'data' / 'delfi' / 'hg19.hic.chr1.6Mb.txt'
    twobit = request.path.parent / 'data' / 'delfi' / 'hg19.chr1.10Mb.2bit'
    twobit = str(twobit)
    blacklist = request.path.parent / 'data' / 'delfi' / 'hg19_darkregion.bed'
    gaps = GenomeGaps.ucsc_hg19()
    
    results = delfi(frag_file, autosomes, bins_file, twobit, blacklist, gaps)


def test_overall_fasta_matches_2bit(request):
    """FASTA-backed DELFI should match the existing 2-bit-backed result."""
    frag_file = request.path.parent / 'data' / 'delfi' / 'hg19.chr1.6Mb.bam'
    autosomes = request.path.parent / 'data' / 'delfi' / 'human.hg19.chr1.6Mb.genome'
    bins_file = request.path.parent / 'data' / 'delfi' / 'hg19.hic.chr1.6Mb.txt'
    twobit = request.path.parent / 'data' / 'delfi' / 'hg19.chr1.10Mb.2bit'
    fasta = request.path.parent / 'data' / 'delfi' / 'hg19.chr1.10Mb.fa'
    blacklist = request.path.parent / 'data' / 'delfi' / 'hg19_darkregion.bed'
    gaps = GenomeGaps.ucsc_hg19()

    twobit_results = delfi(frag_file, autosomes, bins_file, str(twobit), blacklist, gaps)
    fasta_results = delfi(frag_file, autosomes, bins_file, str(fasta), blacklist, gaps)

    pd.testing.assert_frame_equal(twobit_results, fasta_results)


# The two tests below were contributed by D.H.K. (Duco) Gaillard (@DucoG) in
# upstream PR epifluidlab/FinaleToolkit#172, alongside the delfi worker-pool
# speedup they guard.

def _delfi_inputs(request):
    """Shared input paths for the chr1 6Mb hg19 fixture."""
    base = request.path.parent / 'data' / 'delfi'
    return dict(
        input_file=base / 'hg19.chr1.6Mb.bam',
        chrom_sizes=base / 'human.hg19.chr1.6Mb.genome',
        bins_file=base / 'hg19.hic.chr1.6Mb.txt',
        reference_file=str(base / 'hg19.chr1.10Mb.2bit'),
        blacklist_file=base / 'hg19_darkregion.bed',
        gap_file=GenomeGaps.ucsc_hg19(),
    )


def _bam_to_fragfile(bam_path, out_path, bed6=False):
    """Write a sorted, tabix-indexed FinaleDB fragment file from a BAM.

    Fragments are extracted with the BAM reader (no mapq cutoff) so the
    resulting tabix file holds the same fragment coordinates DELFI reads from
    the BAM. ``bed6`` writes a 6-column BED file (with a placeholder name
    field) instead of the 5-column FinaleDB layout. Returns the path to the
    bgzipped, tabix-indexed file.
    """
    records = [
        (contig, int(start), int(stop), int(mapq), '+' if is_forward else '-')
        for contig, start, stop, mapq, is_forward in frag_generator(
            str(bam_path), contig=None, quality_threshold=0)
    ]
    records.sort(key=lambda r: (r[0], r[1], r[2]))
    with open(out_path, 'w') as fh:
        for contig, start, stop, mapq, strand in records:
            if bed6:
                fh.write(f'{contig}\t{start}\t{stop}\t.\t{mapq}\t{strand}\n')
            else:
                fh.write(f'{contig}\t{start}\t{stop}\t{mapq}\t{strand}\n')
    return pysam.tabix_index(str(out_path), preset='bed', force=True)


def test_workers_equivalence(request):
    """delfi(workers=1) must match delfi(workers>1) bit-for-bit.

    Guards against races or parallel-only state in the shared
    worker-process handles set up by the Pool initializer.
    """
    common = dict(
        _delfi_inputs(request),
        no_gc_correct=True, remove_nocov=False, merge_bins=False,
    )
    serial = delfi(**common, workers=1).reset_index(drop=True)
    parallel = delfi(**common, workers=4).reset_index(drop=True)
    pd.testing.assert_frame_equal(serial, parallel)


def test_fragfile_input(request, tmp_path):
    """Tabix-indexed fragment input (.frag.gz and .bed.gz) is supported.

    Regression test for fragment-file input: the worker pool reads tabix
    input via AlignmentWrapper, so it must run without hanging, return
    non-empty output, be invariant to the worker count, and produce identical
    results for the 5-column FinaleDB and 6-column BED layouts.

    Counts are only checked to be *close* to the BAM, not identical: the BAM
    reader fetches each window by read-alignment position while the tabix
    reader fetches by fragment span, so a few window-boundary fragments are
    assigned to adjacent windows. This discrepancy is independent of the
    parallelisation and predates this change.
    """
    inputs = _delfi_inputs(request)
    bam = inputs['input_file']
    shared = {k: v for k, v in inputs.items() if k != 'input_file'}
    common = dict(
        shared, no_gc_correct=True, remove_nocov=False, merge_bins=False,
    )

    frag_gz = _bam_to_fragfile(bam, tmp_path / 'frags.frag', bed6=False)
    bed_gz = _bam_to_fragfile(bam, tmp_path / 'frags.bed', bed6=True)

    frag_serial = delfi(input_file=frag_gz, **common, workers=1).reset_index(drop=True)
    frag_parallel = delfi(input_file=frag_gz, **common, workers=4).reset_index(drop=True)
    bed_parallel = delfi(input_file=bed_gz, **common, workers=4).reset_index(drop=True)

    # tabix input must not hang or silently drop every fragment
    assert frag_parallel.shape[0] > 0
    assert frag_parallel['num_frags'].sum() > 0

    # worker-count invariance for tabix input
    pd.testing.assert_frame_equal(frag_serial, frag_parallel)

    # the two tabix layouts hold identical fragments -> identical output
    pd.testing.assert_frame_equal(frag_parallel, bed_parallel)

    # fragment counts track the BAM closely (see docstring on boundaries)
    bam_result = delfi(input_file=bam, **common, workers=4).reset_index(drop=True)
    bam_total = bam_result['num_frags'].sum()
    frag_total = frag_parallel['num_frags'].sum()
    assert abs(frag_total - bam_total) / bam_total < 0.01


# --------------------------------------------------------------------------- #
# Merge window size (--merge-size / window_size) and bin-size validation
# --------------------------------------------------------------------------- #

def _bins(rows):
    """Build a contig/start/stop DataFrame from ``(contig, start, stop)``."""
    return pd.DataFrame(rows, columns=['contig', 'start', 'stop'])


def _tiled(contig, n, size, inclusive=False, offset=0):
    """``n`` adjacent bins of ``size`` bp, half-open or end-inclusive."""
    return [
        (contig, offset + i * size,
         offset + (i + 1) * size - (1 if inclusive else 0))
        for i in range(n)
    ]


@pytest.mark.parametrize('inclusive', [False, True])
def test_resolve_bin_size_conventions(inclusive):
    """Half-open and end-inclusive (Cristiano HiC style) bins are accepted."""
    bins = _bins(_tiled('chr1', 5, 100, inclusive) + _tiled('chr2', 3, 100, inclusive))
    layout = _resolve_bin_size(bins, 500)
    assert layout.bin_size == 100
    assert layout.inclusive == inclusive
    assert layout.keep.all()


def test_resolve_bin_size_canonical_bins_file(request):
    """The DELFI HiC bins fixture is 100kb end-inclusive."""
    bins_file = request.path.parent / 'data' / 'delfi' / 'hg19.hic.chr1.6Mb.txt'
    bins = pd.read_csv(bins_file, sep='\t', comment='#', usecols=[0, 1, 2],
                       names=['contig', 'start', 'stop'])
    layout = _resolve_bin_size(bins, 5_000_000)
    assert (layout.bin_size, layout.inclusive) == (100_000, True)


@pytest.mark.parametrize('inclusive', [False, True])
def test_resolve_bin_size_drops_short_last_bin(inclusive):
    """A short final bin per contig is flagged for removal, not an error."""
    rows = _tiled('chr1', 4, 100, inclusive) + [('chr1', 400, 450)]
    rows += _tiled('chr2', 2, 100, inclusive) + [('chr2', 200, 230)]
    layout = _resolve_bin_size(_bins(rows), 200)
    assert layout.bin_size == 100
    assert layout.keep.tolist() == [True] * 4 + [False] + [True] * 2 + [False]


@pytest.mark.parametrize('rows', [
    # a short bin in the middle of a contig
    _tiled('chr1', 2, 100) + [('chr1', 200, 250)] + _tiled('chr1', 2, 100, offset=250),
    # a last bin longer than the rest
    _tiled('chr1', 3, 100) + [('chr1', 300, 450)],
    # every contig has a single bin, of different sizes
    [('chr1', 0, 100), ('chr2', 0, 200)],
], ids=['short-middle', 'long-last', 'singletons'])
def test_resolve_bin_size_non_uniform(rows):
    with pytest.raises(ValueError, match='same size'):
        _resolve_bin_size(_bins(rows), 1000)


def test_resolve_bin_size_mixed_convention():
    # Same widths, but one gap of 0 (half-open) and one of 1 (inclusive).
    rows = [('chr1', 0, 100), ('chr1', 100, 200), ('chr1', 201, 301)]
    with pytest.raises(ValueError, match='mix half-open'):
        _resolve_bin_size(_bins(rows), 300)


@pytest.mark.parametrize('window_size', [250, 50, 0, -100])
def test_resolve_bin_size_window_not_multiple(window_size):
    with pytest.raises(ValueError, match='window'):
        _resolve_bin_size(_bins(_tiled('chr1', 5, 100)), window_size)


def test_resolve_bin_size_empty():
    with pytest.raises(ValueError, match='no bins'):
        _resolve_bin_size(_bins([]), 100)


@pytest.fixture
def delfi_100kb_bins(request):
    path = request.path.parent / 'data' / 'delfi' / 'test_delfi_100kb.csv'
    return pd.read_csv(path, dtype={'contig': str, 'start': int, 'stop': int})


def test_merge_bins_window_size(delfi_100kb_bins):
    """window_size controls how many bins are merged per window."""
    merged = delfi_merge_bins(delfi_100kb_bins, window_size=1_000_000)
    for arm, n_merged in (('1p', 1182 // 10), ('1q', 1006 // 10)):
        arm_bins = delfi_100kb_bins[delfi_100kb_bins['arm'] == arm]
        arm_merged = merged[merged['arm'] == arm]
        assert arm_merged.shape[0] == n_merged
        # p-arms are merged from the first bin, q-arms from the last.
        if arm == '1p':
            first = arm_bins.iloc[:10]
            assert arm_merged['start'].iloc[0] == first['start'].iloc[0]
            assert arm_merged['stop'].iloc[0] == first['stop'].iloc[-1]
            assert arm_merged['short'].iloc[0] == first['short'].sum()
        else:
            last = arm_bins.iloc[-10:]
            assert arm_merged['start'].iloc[-1] == last['start'].iloc[0]
            assert arm_merged['stop'].iloc[-1] == last['stop'].iloc[-1]
            assert arm_merged['short'].iloc[-1] == last['short'].sum()


def test_merge_bins_explicit_bin_size_matches_inferred(delfi_100kb_bins):
    inferred = delfi_merge_bins(delfi_100kb_bins, window_size=2_000_000)
    explicit = delfi_merge_bins(
        delfi_100kb_bins, window_size=2_000_000, bin_size=100_000)
    pd.testing.assert_frame_equal(inferred, explicit)


def test_merge_bins_one_bin_per_window(delfi_100kb_bins):
    """window_size == bin size keeps every bin (incl. the first q-arm bin)."""
    merged = delfi_merge_bins(delfi_100kb_bins, window_size=100_000)
    assert merged.shape[0] == delfi_100kb_bins.shape[0]
    assert (merged['start'].to_numpy() == delfi_100kb_bins['start'].to_numpy()).all()


@pytest.mark.parametrize('kwargs', [
    dict(window_size=150_000),
    dict(window_size=150_000, bin_size=100_000),
    dict(window_size=1_000_000, bin_size=0),
])
def test_merge_bins_invalid_window(delfi_100kb_bins, kwargs):
    with pytest.raises(ValueError, match='multiple of the bin size'):
        delfi_merge_bins(delfi_100kb_bins, **kwargs)


def _write_bins(rows, path):
    with open(path, 'w') as fh:
        for contig, start, stop in rows:
            fh.write(f'{contig}\t{start}\t{stop}\n')
    return path


def test_delfi_window_size_reaches_merge(request):
    """delfi(window_size=...) merges the same bins as delfi_merge_bins."""
    common = dict(_delfi_inputs(request), remove_nocov=False)
    unmerged = delfi(**common, merge_bins=False)
    expected = delfi_merge_bins(
        unmerged, window_size=1_000_000, bin_size=100_000)
    merged = delfi(**common, window_size=1_000_000)
    assert merged.shape[0] > delfi(**common).shape[0]
    pd.testing.assert_frame_equal(merged, expected)


def test_delfi_rejects_non_uniform_bins(request, tmp_path):
    bins = _write_bins(
        _tiled('1', 3, 100_000) + [('1', 300_000, 350_000)]
        + _tiled('1', 3, 100_000, offset=350_000),
        tmp_path / 'bins.bed')
    common = dict(_delfi_inputs(request), bins_file=bins)
    with pytest.raises(ValueError, match='same size'):
        delfi(**common)
    # Uniformity is only required when merging.
    delfi(**common, merge_bins=False, no_gc_correct=True, remove_nocov=False)


def test_delfi_rejects_window_not_multiple(request):
    with pytest.raises(ValueError, match='multiple of the bin size'):
        delfi(**_delfi_inputs(request), window_size=1_250_000)


def _cli_delfi_args(request, tmp_path):
    # Fragment-file input: the CLI's BAM-vs-reference contig check would
    # reject the fixture BAM, whose header lists contigs the chr1-only
    # reference lacks.
    inputs = _delfi_inputs(request)
    frag_gz = _bam_to_fragfile(inputs['input_file'], tmp_path / 'frags.frag')
    return [
        'delfi', str(frag_gz), str(inputs['chrom_sizes']),
        inputs['reference_file'], str(inputs['bins_file']),
        '-b', str(inputs['blacklist_file']), '-g', 'hg19',
        '--no-remove-nocov', '-o', str(tmp_path / 'delfi.tsv'),
    ]


def test_cli_merge_size(request, tmp_path):
    args = _cli_delfi_args(request, tmp_path)
    result = CliRunner().invoke(main_cli, args + ['--merge-size', '1000000'])
    assert result.exit_code == 0, result.output
    assert pd.read_csv(tmp_path / 'delfi.tsv', sep='\t').shape[0] == 5


@pytest.mark.parametrize('size', ['0', '-5'])
def test_cli_merge_size_must_be_positive(request, tmp_path, size):
    result = CliRunner().invoke(
        main_cli,
        _cli_delfi_args(request, tmp_path) + ['--merge-size', size])
    assert result.exit_code == 2
    assert '--merge-size' in result.output


def test_merge_bins_whole_contig(delfi_100kb_bins):
    """Without gap annotations (arm == contig), contigs merge 5'->3'."""
    bins = delfi_100kb_bins.assign(contig='chr1', arm='chr1')
    merged = delfi_merge_bins(bins, window_size=1_000_000)
    assert merged.shape[0] == bins.shape[0] // 10
    assert (merged['contig'] == 'chr1').all()
    assert (merged['arm'] == 'chr1').all()
    assert merged['start'].iloc[0] == bins['start'].iloc[0]
    assert merged['stop'].iloc[0] == bins['stop'].iloc[9]
    assert merged['short'].iloc[0] == bins['short'].iloc[:10].sum()


def test_delfi_merges_without_gap_file(request):
    """Merging used to silently return no rows when gap_file was None."""
    common = dict(_delfi_inputs(request), gap_file=None, remove_nocov=False)
    unmerged = delfi(**common, merge_bins=False)
    merged = delfi(**common, window_size=1_000_000)
    expected = delfi_merge_bins(
        unmerged, window_size=1_000_000, bin_size=100_000)
    assert merged.shape[0] == unmerged.shape[0] // 10 > 0
    assert (merged['contig'] == 'chr1').all()
    pd.testing.assert_frame_equal(merged, expected)


# --------------------------------------------------------------------------- #
# gap_file resolution (-g/--gap-file)
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize('name,factory', [
    ('hg19', GenomeGaps.ucsc_hg19),
    ('b37', GenomeGaps.b37),
    ('human_g1k_v37', GenomeGaps.b37),
    ('hg38', GenomeGaps.hg38),
    ('GRCh38', GenomeGaps.hg38),
])
def test_resolve_gaps_genome_name(name, factory):
    np.testing.assert_array_equal(_resolve_gaps(name).gaps, factory().gaps)


@pytest.mark.parametrize('as_path', [False, True])
def test_resolve_gaps_bed_file(tmp_path, as_path):
    bed = tmp_path / 'gaps.bed'
    GenomeGaps.ucsc_hg19().to_bed(bed)
    gaps = _resolve_gaps(bed if as_path else str(bed))
    np.testing.assert_array_equal(gaps.gaps, GenomeGaps.ucsc_hg19().gaps)


def test_resolve_gaps_existing_file_wins_over_name(tmp_path, monkeypatch):
    """A file literally named like a genome is read as a file."""
    monkeypatch.chdir(tmp_path)
    GenomeGaps.b37().to_bed('hg19')
    np.testing.assert_array_equal(_resolve_gaps('hg19').gaps, GenomeGaps.b37().gaps)


def test_resolve_gaps_unknown():
    with pytest.raises(FileNotFoundError, match='known genome name'):
        _resolve_gaps('not-a-genome')
