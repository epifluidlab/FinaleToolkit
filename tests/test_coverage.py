"""
Tests for finaletoolkit.frag.coverage module and associated
subcommands, which calculate fragment coverage.
"""

import os
import filecmp
import difflib

import pytest

from finaletoolkit.frag import coverage, single_coverage

class TestSingleCoverage:
    def test_coverage(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        chrom, start, stop, name, cov = single_coverage(bam, '12', 0, None,
                                                        quality_threshold=0)

        assert chrom == '12'
        assert start == 0
        assert cov == pytest.approx(17)

    def test_coverage_interval(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        chrom, start, stop, name, cov = single_coverage(
            bam, '12', 34443000, 34447000, quality_threshold=0)

        assert chrom == '12'
        assert start == 34443000
        assert stop == 34447000
        assert cov == pytest.approx(17)

    def test_coverage_interval_midpoints(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        chrom, start, stop, name, cov = single_coverage(
            bam, '12', 34443400, 34443600, quality_threshold=0)

        assert chrom == '12'
        assert start == 34443400
        assert stop == 34443600
        assert cov == pytest.approx(2)

class TestCoverage:
    def test_coverage_normalize(self, request):
        input_file = request.path.parent / 'data' / '12.3444.b37.frag.gz'
        intervals = request.path.parent / 'data' / 'intervals.bed'
        results = coverage(
            input_file,
            intervals,
            "-",
            scale_factor=1.,
            normalize=True,
            verbose=True)
        for i in range(2):
            chrom, start, stop, name, cov = results[i]
            if i == 0:  
                assert chrom == '12'
                assert start == 34443118
                assert stop == 34443538
                assert name == '.'
                assert cov == pytest.approx(4/16)
            elif i == 1:
                assert chrom == '12'
                assert start == 34444968
                assert stop == 34446115
                assert name == '.'
                assert cov == pytest.approx(7/16)	

    def test_coverage_no_normalize(self, request):
        input_file = request.path.parent / 'data' / '12.3444.b37.frag.gz'
        intervals = request.path.parent / 'data' / 'intervals.bed'
        results = coverage(input_file, intervals, "-", normalize=False,
                           intersect_policy='midpoint',
                           scale_factor=1.)
        for i in range(2):
            chrom, start, stop, name, cov = results[i]
            if i == 0:  
                assert chrom == '12'
                assert start == 34443118
                assert stop == 34443538
                assert name == '.'
                assert cov == pytest.approx(4)
            elif i == 1:
                assert chrom == '12'
                assert start == 34444968
                assert stop == 34446115
                assert name == '.'
                assert cov == pytest.approx(7)	
