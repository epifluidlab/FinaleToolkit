"""
Tests for finaletoolkit.frag.end_motifs module, which calculates
end-motifs and motif diversity scores genomewide and over intervals.
"""
import pytest
import os
import filecmp

from finaletoolkit.frag.end_motifs import *

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
    
    def test_tsv(self, request, tmp_path):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        ref_file = request.path.parent / 'data' / 'b37.chr12.2bit'
        motifs = end_motifs(
            bam, ref_file, both_strands=True, quality_threshold=0)
        
        dest = tmp_path / "results.tsv"
        dif_file = request.path.parent / 'data' / 'end_motifs_dif.tsv'
        motifs.to_tsv(dest)
        assert filecmp.cmp(dest, dif_file)

        
class TestInvervalEndMotifs:
    def test_end_motifs(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        ref_file = request.path.parent / 'data' / 'b37.chr12.2bit'
        motifs = interval_end_motifs(
            bam, ref_file, [('12', 34440000, 34450000)], both_strands=True,
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
            bam, ref_file, [('12', 34442500, 34446500)], both_strands=True,
            quality_threshold=0)
        
        assert motifs.motif_diversity_score()[0][1] == pytest.approx(
            0.5844622669209985, 0.1)


class TestCLIEntryPoint:
    """
    Test all CLI subcommands related to end_motifs and MDS, genomewide
    and intervals.
    """
    def test_end_motif(self, request):
        exit_status = os.system('finaletoolkit end-motifs --help')
        assert exit_status == 0

    def test_mds(self, request):
        exit_status = os.system('finaletoolkit mds --help')
        assert exit_status == 0

    def test_interval_end_motif(self, request):
        exit_status = os.system('finaletoolkit interval-end-motifs --help')
        assert exit_status == 0

    def test_interval_mds(self, request):
        exit_status = os.system('finaletoolkit interval-mds --help')
        assert exit_status == 0
    