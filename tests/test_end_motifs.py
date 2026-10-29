"""
Tests for finaletoolkit.frag.end_motifs module, which calculates
end-motifs and motif diversity scores genomewide and over intervals.
"""

import os
import filecmp
import difflib

import numpy as np
import pytest
from click.testing import CliRunner

from finaletoolkit.cli.main_cli import main_cli
from finaletoolkit.frag import end_motifs, EndMotifFreqs, EndMotifsIntervals, interval_end_motifs
from finaletoolkit.frag import region_breakpoint_motifs, region_end_motifs
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

class TestFragmentLengthFilter:
    """min/max length bounds reach the per-region counters (issue: they were
    silently ignored after the AlignmentWrapper refactor in v0.12.0).

    ``12.3444.b37.bam`` holds 17 fragments, 121-205 bp; 8 are 150-170 bp and
    two are exactly 161 bp. With both strands, each fragment adds two ends.
    """

    REGION = ('12', 34_000_000, 35_000_000)

    @staticmethod
    def _files(request):
        data = request.path.parent / 'data'
        return str(data / '12.3444.b37.bam'), str(data / 'b37.chr12.2bit')

    @pytest.mark.parametrize('low, high, n_frags', [
        (None, None, 17),
        (150, 170, 8),
        (161, 161, 2),  # bounds are inclusive
        (300, None, 0),
    ])
    def test_region_end_motifs(self, request, low, high, n_frags):
        bam, ref = self._files(request)
        counts = region_end_motifs(
            bam, *self.REGION, ref, fraction_low=low, fraction_high=high,
            quality_threshold=0)
        assert sum(counts.values()) == 2 * n_frags

    @pytest.mark.parametrize('low, high, n_frags', [
        (10, 600, 17),
        (150, 170, 8),
        (300, None, 0),
    ])
    def test_region_breakpoint_motifs(self, request, low, high, n_frags):
        bam, ref = self._files(request)
        counts = region_breakpoint_motifs(
            bam, *self.REGION, ref, fraction_low=low, fraction_high=high,
            quality_threshold=0)
        assert sum(counts.values()) == 2 * n_frags

    def test_interval_end_motifs(self, request):
        bam, ref = self._files(request)
        motifs = interval_end_motifs(
            bam, ref, [(*self.REGION, 'all')], min_length=150,
            max_length=170, quality_threshold=0)
        (_, counts), = motifs.intervals
        assert sum(counts.values()) == 16

    def test_genomewide_end_motifs(self, request):
        bam, ref = self._files(request)
        motifs = end_motifs(
            bam, ref, min_length=161, max_length=161, quality_threshold=0)
        # Two fragments -> four ends, so every frequency is a multiple of 1/4.
        freqs = np.array(motifs.frequencies())
        assert freqs.sum() == pytest.approx(1)
        assert freqs * 4 == pytest.approx(np.round(freqs * 4))

    def test_cli_max_length(self, request, tmp_path):
        # Fragment file: the BAM header has contigs missing from the chr12-only
        # reference, which the CLI rejects. Same 17 fragments.
        _, ref = self._files(request)
        frags = str(request.path.parent / 'data' / '12.3444.b37.frag.gz')
        intervals = tmp_path / 'all.bed'
        intervals.write_text('\t'.join(map(str, self.REGION)) + '\tall\n')
        out = tmp_path / 'out.tsv'
        result = CliRunner().invoke(main_cli, [
            'interval-end-motifs', frags, ref, str(intervals), '-q', '0',
            '--max-length', '170', '-o', str(out)])
        assert result.exit_code == 0, result.output
        motifs = EndMotifsIntervals.from_file(str(out), 0, '\t', 0)
        # 150-170 bp plus 121, 135, 137, 147 bp -> 12 fragments.
        assert motifs.total_counts == [24]


class TestMidpointMembership:
    """A fragment belongs to a region only if its midpoint lies in
    ``[start, stop)``, as before the v0.12.0 AlignmentWrapper refactor.

    Five fragments overlap 12:34443118-34443538; the last (34443483-34443660,
    midpoint 34443571) has its midpoint past the end. A split at 34443200
    cuts through the first two fragments (and their read1 alignments, which
    is what a BAM index query matches on).
    """

    @staticmethod
    def _files(request):
        data = request.path.parent / 'data'
        return str(data / '12.3444.b37.bam'), str(data / 'b37.chr12.2bit')

    def test_overlap_without_midpoint_is_excluded(self, request):
        bam, ref = self._files(request)
        counts = region_end_motifs(
            bam, '12', 34443118, 34443538, ref, quality_threshold=0)
        assert sum(counts.values()) == 2 * 4

    @pytest.mark.parametrize('input_name', [
        '12.3444.b37.bam', '12.3444.b37.frag.gz',
    ])
    @pytest.mark.parametrize('region_motifs', [
        region_end_motifs, region_breakpoint_motifs,
    ])
    def test_adjacent_regions_partition_fragments(
        self, request, region_motifs, input_name
    ):
        _, ref = self._files(request)
        frags = str(request.path.parent / 'data' / input_name)

        def total(start, stop):
            return sum(region_motifs(
                frags, '12', start, stop, ref, quality_threshold=0).values())

        whole = total(34443118, 34443700)
        assert whole == 2 * 5
        assert total(34443118, 34443200) + total(34443200, 34443700) == whole
