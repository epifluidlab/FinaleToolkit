"""
Tests for finaletoolkit.frag.wps, finaletoolkit.frag.adjust_wps, and
finaletoolkit.frag.multi_wps
"""

import os
import filecmp
import difflib

import pytest
import numpy as np

from finaletoolkit.frag import wps, multi_wps

class TestWPS:
    def test_lwps(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        results = wps(bam, '12', 34444145, 34444155, quality_threshold=0)

        assert np.all(results['contig'] == '12')
        assert np.all(
            results['start'] == np.arange(34444145, 34444155)), str(results)
        assert np.all(
            results['wps'] == [-1,-1,-1,-1,-1,1,1,1,1,1]), str(results)