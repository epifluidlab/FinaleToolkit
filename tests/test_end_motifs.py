"""
Tests for finaletoolkit.frag.end_motifs module, which calculates
end-motifs and motif diversity scores genomewide and over intervals.
"""

import os
import filecmp
import difflib

import numpy as np
import pytest

from finaletoolkit.frag import end_motifs, EndMotifFreqs, EndMotifsIntervals, interval_end_motifs
from finaletoolkit.utils import gen_kmers

class TestGenomewideEndMotifs:
    def test_end_motifs(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        ref_file = request.path.parent / 'data' / 'b37.chr12.2bit'
        motifs = end_motifs(
            bam, ref_file, both_strands=True, quality_threshold=0)
        expected_freqs = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.058823529411764705,
                          0.058823529411764705,0.058823529411764705,0.0,0.0,
                          0.0,0.0,0.0,0.058823529411764705,
                          0.058823529411764705,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.029411764705882353,0.0,
                          0.0,0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.058823529411764705,0.0,0.0,0.029411764705882353,
                          0.0,0.0,0.0,0.029411764705882353,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.029411764705882353,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.058823529411764705,
                          0.0,0.029411764705882353,0.0,0.029411764705882353,
                          0.0,0.0,0.0,0.0,0.0,0.029411764705882353,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        for freq, expected in zip(motifs.frequencies(), expected_freqs):
            assert freq == pytest.approx(expected, 0.1)

    def test_mds(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        ref_file = request.path.parent / 'data' / 'b37.chr12.2bit'
        motifs = end_motifs(
            bam, ref_file, both_strands=True, quality_threshold=0)
        
        assert motifs.motif_diversity_score() == pytest.approx(
            0.5844622669209985, 0.1)
    
    def test_to_tsv(self, request, tmp_path):
        """
        Dif output with expected output file.
        """
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        ref_file = request.path.parent / 'data' / 'b37.chr12.2bit'
        motifs = end_motifs(
            bam, ref_file, both_strands=True, quality_threshold=0)
        
        dest = tmp_path / "results.tsv"
        dif_file = request.path.parent / 'data' / 'end_motifs' / 'end_motifs_dif.tsv'
        motifs.to_tsv(dest)
        filecmp.clear_cache()
        assert filecmp.cmp(dest, dif_file)

    def test_from_file(self, request, tmp_path):
        """
        Compare EndMotifFreq with EndMotifFreq created from 
        """
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        ref_file = request.path.parent / 'data' / 'b37.chr12.2bit'
        motifs = end_motifs(
            bam, ref_file, both_strands=True, quality_threshold=0)
        
        dest = tmp_path / "results.tsv"
        motifs.to_tsv(dest)
        tsv_motifs = EndMotifFreqs.from_file(dest, 0)

        assert motifs.freq_dict == pytest.approx(tsv_motifs.freq_dict) 


class TestInvervalEndMotifs:
    def test_end_motifs(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        ref_file = request.path.parent / 'data' / 'b37.chr12.2bit'
        motifs = interval_end_motifs(
            bam, ref_file, [('12', 34440000, 34450000, '.')], both_strands=True,
            quality_threshold=0)
        expected_freqs = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.058823529411764705,
                          0.058823529411764705,0.058823529411764705,0.0,0.0,
                          0.0,0.0,0.0,0.058823529411764705,
                          0.058823529411764705,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.029411764705882353,0.0,
                          0.0,0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.058823529411764705,0.0,0.0,0.029411764705882353,
                          0.0,0.0,0.0,0.029411764705882353,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.0,0.029411764705882353,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.058823529411764705,
                          0.0,0.029411764705882353,0.0,0.029411764705882353,
                          0.0,0.0,0.0,0.0,0.0,0.029411764705882353,0.0,0.0,
                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                          0.029411764705882353,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        for freq, expected in zip(motifs.intervals[0][1].values(), expected_freqs):
            assert freq == pytest.approx(expected/0.029411764705882353, 0.1)

    def test_mds(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        ref_file = request.path.parent / 'data' / 'b37.chr12.2bit'
        motifs = interval_end_motifs(
            bam, ref_file, [('12', 34442500, 34446500, '.')], both_strands=True,
            quality_threshold=0)
        
        assert motifs.motif_diversity_score()[0][1] == pytest.approx(
            0.5844622669209985, 0.1)
        
    def test_to_tsv(self, request, tmp_path):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        ref_file = request.path.parent / 'data' / 'b37.chr12.2bit'
        motifs = interval_end_motifs(
            bam, ref_file, [('12', 34440000, 34450000, '.')], both_strands=True,
            quality_threshold=0)
        
        dest = tmp_path / "interval_results.tsv"
        dif_file = request.path.parent / 'data' / 'end_motifs' / 'end_motifs_intervals_dif.tsv'

        motifs.to_tsv(dest, False)
        differ = difflib.Differ()
        with open(dest) as actual, open(dif_file) as expected:
            difs = differ.compare(actual.readlines(), expected.readlines())
            dif_text = '\n'.join(difs)
        filecmp.clear_cache()
        assert filecmp.cmp(dest, dif_file), (
            f"Generated file is different from expected: \n{dif_text}")

    def test_from_file(self, request, tmp_path):
        pass

    
class TestCLIRegionalMDS:
    def test_regional_mds(self, request, tmp_path):
        src_tsv = request.path.parent / 'data' / 'end_motifs' / 'end_motifs_intervals_dif.tsv'
        dest = tmp_path / "regional_mds.tsv"
        os.system(f'finaletoolkit regional-mds {src_tsv} {dest}')
        with open(dest) as file:
            line = file.readline()
            chrom, start, stop, name, mds = line.split()
        assert float(mds) == pytest.approx(0.5844622669209985, 0.01)

    def test_regional_mds_miller_madow(self, request, tmp_path):
        """--miller-madow raises rMDS by the analytic (m-1)/(2N) term."""
        src_tsv = request.path.parent / 'data' / 'end_motifs' / 'end_motifs_intervals_dif.tsv'
        dest = tmp_path / "regional_mds_mm.tsv"
        os.system(f'finaletoolkit regional-mds {src_tsv} {dest} --miller-madow')
        with open(dest) as file:
            chrom, start, stop, name, mds = file.readline().split()
        # The region holds N=34 fragment ends across m=27 observed 4-mers, so
        # the correction is (27 - 1) / (2 * 34) / log(4**4).
        expected = 0.5844622669209985 + 26 / 68 / np.log(256)
        assert float(mds) == pytest.approx(expected, 0.01)


class TestRegionalMDSMillerMadow:
    """Miller-Madow bias correction on the rMDS API."""

    @staticmethod
    def _intervals(request):
        src_tsv = (
            request.path.parent / 'data' / 'end_motifs'
            / 'end_motifs_intervals_dif.tsv'
        )
        return EndMotifsIntervals.from_file(str(src_tsv), 30, '\t', 0)

    def test_off_by_default(self, request):
        motifs = self._intervals(request)
        assert (
            motifs.motif_diversity_score()
            == motifs.motif_diversity_score(miller_madow=False)
        )

    def test_correction_matches_formula(self, request):
        motifs = self._intervals(request)
        (_, plain), = motifs.motif_diversity_score()
        (_, corrected), = motifs.motif_diversity_score(miller_madow=True)

        _, kmers = motifs.intervals[0]
        counts = np.array(list(kmers.values()))
        n = counts.sum()
        occupied = int((counts > 0).sum())

        assert corrected > plain
        assert corrected - plain == pytest.approx(
            (occupied - 1) / (2 * n) / np.log(4 ** motifs.k)
        )

    def test_frequency_table_round_trip(self, request, tmp_path):
        """N survives a to_tsv/from_file round trip through frequencies.

        ``to_tsv(calc_freq=True)`` normalizes the per-motif values away, so the
        correction has to read N from the retained ``count`` column rather than
        summing the row (which would give 1).
        """
        motifs = self._intervals(request)
        dest = tmp_path / "freqs.tsv"
        motifs.to_tsv(dest, calc_freq=True)
        reloaded = EndMotifsIntervals.from_file(str(dest), 30, '\t', 0)

        (_, from_counts), = motifs.motif_diversity_score(miller_madow=True)
        (_, from_freqs), = reloaded.motif_diversity_score(miller_madow=True)
        assert from_freqs == pytest.approx(from_counts)

    def test_empty_region_is_nan(self):
        """A region with no fragments yields NaN rather than dividing by zero."""
        kmers = dict.fromkeys(gen_kmers(4, "ACGT"), 0.0)
        motifs = EndMotifsIntervals([(("12", 0, 100, "."), kmers)], 4, 30)
        (_, score), = motifs.motif_diversity_score(miller_madow=True)
        assert np.isnan(score)