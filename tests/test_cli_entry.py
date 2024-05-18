"""
Tests for main_cli and entry points.
"""
import os

import pytest

class TestCLIEntryPoint:
    """
    Test all CLI subcommands related to end_motifs and MDS, genomewide
    and intervals.
    """
    def test_coverage(self, request):
        exit_status = os.system('finaletoolkit coverage --help')
        assert exit_status == 0

    def test_frag_length_bins(self, request):
        exit_status = os.system('finaletoolkit frag-length-bins --help')
        assert exit_status == 0

    def test_frag_length_intervals(self, request):
        exit_status = os.system('finaletoolkit frag-length-intervals --help')
        assert exit_status == 0

    def test_cleavage_profile(self, request):
        exit_status = os.system('finaletoolkit cleavage-profile --help')
        assert exit_status == 0

    def test_wps(self, request):
        exit_status = os.system('finaletoolkit wps --help')
        assert exit_status == 0

    def test_agg_wps(self, request):
        exit_status = os.system('finaletoolkit adjust-wps --help')
        assert exit_status == 0

    def test_delfi(self, request):
        exit_status = os.system('finaletoolkit delfi --help')
        assert exit_status == 0

    def test_delfi_gc_correct(self, request):
        exit_status = os.system('finaletoolkit delfi-gc-correct --help')
        assert exit_status == 0

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

    def test_filter_bam(self, request):
        exit_status = os.system('finaletoolkit filter-bam --help')
        assert exit_status == 0

    def test_agg_bw(self, request):
        exit_status = os.system('finaletoolkit agg-bw --help')
        assert exit_status == 0

    def test_gap_bed(self, request):
        exit_status = os.system('finaletoolkit gap-bed --help')
        assert exit_status == 0
    