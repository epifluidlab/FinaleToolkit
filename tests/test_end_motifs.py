"""
Tests for finaletoolkit.frag.end_motifs module, which calculates
end-motifs and motif diversity scores genomewide and over intervals.
"""

import os
import filecmp
import difflib

import pytest

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
    
    def test_to_tsv(self, request, tmp_path):
        """
        Dif output with expected output file.
        """
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        ref_file = request.path.parent / 'data' / 'b37.chr12.2bit'
        motifs = end_motifs(
            bam, ref_file, both_strands=True, quality_threshold=0)
        
        dest = tmp_path / "results.tsv"
        dif_file = request.path.parent / 'data' / 'end_motifs_dif.tsv'
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
        dif_file = request.path.parent / 'data' / 'end_motifs_intervals_dif.tsv'

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
    
class TestCLIIntervalMDS:
    def test_interval_mds(self, request, tmp_path):
        src_tsv = request.path.parent / 'data' / 'end_motifs_intervals_dif.tsv'
        dest = tmp_path / "interval_mds.tsv"
        os.system(f'finaletoolkit interval-mds {src_tsv} {dest}')
        with open(dest) as file:
            line = file.readline()
            chrom, start, stop, name, mds = line.split() 
        assert float(mds) == pytest.approx(0.5844622669209985, 0.01)