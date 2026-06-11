"""
Tests for modules related to calculating DELFI, including 
finaletoolkit.frag.delfi_merge_bins, finaletoolkit.frag.delfi, and
finaletoolkit.frag.delfi_gc_correct.
"""

import pytest

import pandas as pd
import pysam

from finaletoolkit.frag import delfi_merge_bins
from finaletoolkit.frag import delfi
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

    Fragments are extracted with the BAM reader (no mapq cutoff) so that the
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

    Regression test for fragment-file input being dropped: the worker pool
    reads tabix input via AlignmentWrapper (not pysam.AlignmentFile), so it
    must run without hanging, return non-empty output, be invariant to the
    worker count, and produce identical results for the 5-column FinaleDB
    and 6-column BED layouts.

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

