"""
Tests for finaletoolkit.utils.agg_bw().
"""

import pytest

from finaletoolkit.utils._agg_bw import agg_bw

class TestAggBigWig:
    def test_agg_bw(self, request, tmp_path):
        """Test if it runs"""
        bw = request.path.parent / 'data' / 'test.bw'
        intervals = request.path.parent / 'data' / 'bw_test.bed'
        dest = tmp_path / "out.wig"
        scores = agg_bw(bw, intervals, dest, 0)

        assert scores == pytest.approx([0.,0.,0.,0.,0.])

    def test_median(self, request, tmp_path):
        """Test median filter"""
        bw = request.path.parent / 'data' / 'test.bw'
        intervals = request.path.parent / 'data' / 'bw_test.bed'
        dest = tmp_path / "out.wig"
        scores = agg_bw(bw, intervals, dest, 2)

        assert scores == pytest.approx([1.,2.,3.])