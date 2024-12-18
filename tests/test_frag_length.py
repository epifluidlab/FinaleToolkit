"""
Tests for finaletoolkit.frag.frag_length module and associated
subcommands, which calculate fragment coverage.
"""

import os
import filecmp
import difflib

import pytest
import numpy as np

from finaletoolkit.frag import frag_length, frag_length_bins, frag_length_intervals

class TestFragLength:
    def test_frag_lengths(self, request):
        frag_file = request.path.parent / 'data' / '12.3444.b37.frag.gz'
        contig = "12"
        start=34443119
        stop=34443538
        lengths = frag_length(frag_file, contig=contig,start=start,stop=stop)
        expected = [166, 161, 180, 177]
        assert np.any(np.equal(lengths, expected))

class TestFragLengthBins:
    def test_default(self, request):
        frag_file = request.path.parent / 'data' / '12.3444.b37.frag.gz'
        contig = "12"
        start=34443119
        stop=34443538
        bins, counts = frag_length_bins(frag_file, contig=contig,start=start,stop=stop)
        expected_bins = [166, 161, 180, 177]
        for bin in expected_bins:
            assert np.isin(bin, bins)
        for count in counts:
            assert count == 1 or count == 0

class TestFragLengthIntervals:
    def test_default(self, request):
        """
        Placeholder test to see if procedure executes
        """
        frag_file = request.path.parent / 'data' / '12.3444.b37.frag.gz'
        intervals = request.path.parent / 'data' / 'intervals.bed'
        results = frag_length_intervals(frag_file, intervals)