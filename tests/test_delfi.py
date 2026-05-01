"""
Tests for modules related to calculating DELFI, including 
finaletoolkit.frag.delfi_merge_bins, finaletoolkit.frag.delfi, and
finaletoolkit.frag.delfi_gc_correct.
"""

import pytest

import pandas as pd

from finaletoolkit.frag import delfi_merge_bins
from finaletoolkit.frag import delfi
from finaletoolkit.genome.gaps import GenomeGaps


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


def test_workers_equivalence(request):
    """delfi(workers=1) must produce identical output to delfi(workers>1).

    Guards against races or parallel-only state in the worker pool.
    """
    frag_file = request.path.parent / 'data' / 'delfi' / 'hg19.chr1.6Mb.bam'
    autosomes = request.path.parent / 'data' / 'delfi' / 'human.hg19.chr1.6Mb.genome'
    bins_file = request.path.parent / 'data' / 'delfi' / 'hg19.hic.chr1.6Mb.txt'
    twobit = str(request.path.parent / 'data' / 'delfi' / 'hg19.chr1.10Mb.2bit')
    blacklist = request.path.parent / 'data' / 'delfi' / 'hg19_darkregion.bed'
    gaps = GenomeGaps.ucsc_hg19()

    common = dict(
        input_file=frag_file, chrom_sizes=autosomes, bins_file=bins_file,
        reference_file=twobit, blacklist_file=blacklist, gap_file=gaps,
        no_gc_correct=True, remove_nocov=False, merge_bins=False,
    )
    serial = delfi(workers=1, **common)
    parallel = delfi(workers=4, **common)

    # bit-identical on every column
    cols = ['contig', 'start', 'stop', 'arm', 'short', 'long', 'gc', 'num_frags']
    for col in cols:
        s = serial[col].fillna(-1).reset_index(drop=True)
        p = parallel[col].fillna(-1).reset_index(drop=True)
        assert s.equals(p), f"workers=1 vs workers=4 differ on column {col}"
