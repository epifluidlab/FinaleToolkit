"""
Tests for finaletoolkit.frag.coverage module and associated
subcommands, which calculate fragment coverage.
"""

import os
import filecmp
import difflib

import pytest

from finaletoolkit.frag.coverage import *

class TestSingleCoverage:
    def test_coverage(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        chrom, start, stop, name, cov = single_coverage(bam, '12', 0, None, quality_threshold=0)

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
    def coverage(self, request):
        pass