"""
Tests for finaletoolkit.frag.coverage module and associated
subcommands, which calculate fragment coverage.
"""

import os
import filecmp
import difflib

import pytest

from finaletoolkit.frag.coverage import *
from finaletoolkit.frag.coverage import _single_coverage

class TestSingleCoverage:
    def test_coverage(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        chrom, start, stop, name, cov = _single_coverage(bam, '12', 0, None, '.', None, None, 1., "midpoint", 0, False)

        assert chrom == '12'
        assert start == 0
        assert cov == pytest.approx(17)

    def test_coverage_interval(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        chrom, start, stop, name, cov = _single_coverage(
            bam, '12', 34443000, 34447000, '.', None, None, 1., "midpoint", 0, False)

        assert chrom == '12'
        assert start == 34443000
        assert stop == 34447000
        assert cov == pytest.approx(17)

    def test_coverage_interval_midpoints(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        chrom, start, stop, name, cov = _single_coverage(
            bam, '12', 34443400, 34443600, '.', None, None, 1., "midpoint", 0, False)

        assert chrom == '12'
        assert start == 34443400
        assert stop == 34443600
        assert cov == pytest.approx(2)

    def test_deprecated_coverage(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        chrom, start, stop, name, cov = single_coverage(bam, '12', 0, None, '.', None, None, 1., "midpoint", 0, False)

        assert chrom == '12'
        assert start == 0
        assert cov == pytest.approx(17)

class TestCoverage:
    def test_func(self, request):
        input_file = request.path.parent / 'data' / '12.3444.b37.frag.gz'
        intervals = request.path.parent / 'data' / 'intervals.bed'
        results = coverage(input_file, intervals, "-", lazy=False)
        print(len(results))
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
                assert cov == pytest.approx(8)

    def test_lazy_func(self, request):
        input_file = request.path.parent / 'data' / '12.3444.b37.frag.gz'
        intervals = request.path.parent / 'data' / 'intervals.bed'
        results = coverage(input_file, intervals, "-", lazy=True)
        print(len(results))
        for result in results:
            chrom, start, stop, name, cov = result
            assert chrom == '12'
            assert start == 34443118 or start == 34444968
            assert stop == 34443538 or stop == 34446115
            assert name == '.' or name == '.'
            assert cov == pytest.approx(312500.0) or cov == pytest.approx(437500.0)
