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

class TestCoverage:
    def coverage(self, request):
        pass