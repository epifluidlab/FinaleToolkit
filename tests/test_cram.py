"""
Integration tests verifying that CRAM files produce identical results to
BAM files when reference_file is supplied. Tests use the existing DELFI
test data (hg19.chr1.6Mb.bam + hg19.chr1.10Mb.fa).

All tests are skipped automatically when samtools is not on PATH.
"""

from __future__ import annotations
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from finaletoolkit.frag import single_coverage, frag_length_bins, wps, delfi
from finaletoolkit.genome.gaps import GenomeGaps


pytestmark = pytest.mark.skipif(
    shutil.which("samtools") is None,
    reason="samtools not on PATH"
)

DATA = Path(__file__).parent / "data" / "delfi"
BAM = DATA / "hg19.chr1.6Mb.bam"
FASTA = DATA / "hg19.chr1.10Mb.fa"
CRAM = DATA / "hg19.chr1.6Mb.cram"


@pytest.fixture(scope="session", autouse=True)
def cram_file():
    """Create CRAM + index from the DELFI BAM and FASTA once per session."""
    if not CRAM.exists():
        subprocess.run(
            ["samtools", "view", "-C", "-T", str(FASTA), str(BAM), "-o", str(CRAM)],
            check=True,
        )
        subprocess.run(["samtools", "index", str(CRAM)], check=True)
    yield CRAM
    # leave the CRAM in place for reuse across runs


class TestSingleCoverageCRAM:
    def test_cram_matches_bam(self, cram_file):
        bam_result = single_coverage(BAM, "chr1", 0, 5000000, quality_threshold=0)
        cram_result = single_coverage(
            cram_file, "chr1", 0, 5000000,
            quality_threshold=0,
            reference_file=FASTA,
        )
        assert bam_result[4] == cram_result[4], (
            f"Coverage differs: BAM={bam_result[4]}, CRAM={cram_result[4]}"
        )


class TestFragLengthBinsCRAM:
    def test_cram_matches_bam(self, cram_file):
        bam_bins, bam_counts = frag_length_bins(
            BAM, contig="chr1", start=0, stop=5000000, quality_threshold=0
        )
        cram_bins, cram_counts = frag_length_bins(
            cram_file, contig="chr1", start=0, stop=5000000,
            quality_threshold=0,
            reference_file=FASTA,
        )
        np.testing.assert_array_equal(bam_bins, cram_bins)
        np.testing.assert_array_equal(bam_counts, cram_counts)


class TestWpsCRAM:
    def test_cram_matches_bam(self, cram_file):
        chrom_size = 249250621  # hg19 chr1
        bam_scores = wps(
            BAM, "chr1", 1000000, 1001000, chrom_size,
            quality_threshold=0
        )
        cram_scores = wps(
            cram_file, "chr1", 1000000, 1001000, chrom_size,
            quality_threshold=0,
            reference_file=FASTA,
        )
        np.testing.assert_array_equal(bam_scores, cram_scores)


class TestDelfiCRAM:
    def test_cram_matches_bam(self, cram_file):
        import pandas as pd

        autosomes = DATA / "human.hg19.chr1.6Mb.genome"
        bins_file = DATA / "hg19.hic.chr1.6Mb.txt"
        blacklist = DATA / "hg19_darkregion.bed"
        gaps = GenomeGaps.ucsc_hg19()
        fasta = str(FASTA)

        # Both use FASTA as reference_file (4th positional). BAM doesn't
        # need it for reading, but CRAM does — ensuring identical GC inputs.
        bam_result = delfi(BAM, autosomes, bins_file, fasta, blacklist, gaps)
        cram_result = delfi(cram_file, autosomes, bins_file, fasta, blacklist, gaps)

        pd.testing.assert_frame_equal(bam_result, cram_result)
