"""
Tests for finaletoolkit.frag.end_motifs module and associated
subcommands, which calculate cleavage profile.
"""

import os
import filecmp
import difflib

import pytest
import numpy as np

from finaletoolkit.frag import cleavage_profile

class TestIntervalCleavageProfile:
    def test_cpg_34443200(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        results = cleavage_profile(
            bam, 133851895, '12', 34443200, 34443201, 5, 5,
            quality_threshold=0)

        assert np.all(results['contig'] == '12')
        assert np.all(results['pos'] == np.arange(34443195, 34443206)), str(results)
        assert np.all(results['proportion'] == np.zeros(11, np.int64)), str(results)
        
