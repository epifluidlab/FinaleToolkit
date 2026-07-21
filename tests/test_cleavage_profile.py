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
from finaletoolkit.frag._cleavage_profile import _coverage_and_ends

class TestIntervalCleavageProfile:
    def test_cpg_34443200(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        results = cleavage_profile(
            bam, 133851895, '12', 34443200, 34443201, 5, 5,
            quality_threshold=0)

        assert np.all(results['contig'] == '12')
        assert np.all(results['pos'] == np.arange(34443195, 34443206)), str(results)
        assert np.all(results['proportion'] == np.zeros(11, np.int64)), str(results)


def _brute_force_coverage_and_ends(starts, stops, strands, adj_start, adj_stop):
    """Reference implementation: the original (n_fragments, n_positions)
    broadcast-matrix algorithm that _coverage_and_ends replaces."""
    positions = np.arange(adj_start, adj_stop)

    fragwise_overlaps = np.logical_and(
        np.greater_equal(positions[np.newaxis], starts[:, np.newaxis]),
        np.less(positions[np.newaxis], stops[:, np.newaxis]),
    )
    depth = np.sum(fragwise_overlaps, axis=0)

    forward_ends = np.logical_and(
        np.equal(positions[np.newaxis], starts[:, np.newaxis]),
        strands[:, np.newaxis],
    )
    reverse_ends = np.logical_and(
        np.equal(positions[np.newaxis], stops[:, np.newaxis]),
        np.logical_not(strands[:, np.newaxis]),
    )
    ends = np.sum(np.logical_or(forward_ends, reverse_ends), axis=0)

    return depth, ends


class TestCoverageAndEndsEquivalence:
    """_coverage_and_ends (diff-array) must match the brute-force
    broadcast-matrix algorithm it replaces, on both random fragment sets
    and boundary edge cases."""

    @pytest.mark.parametrize("seed", range(25))
    def test_random_fragments(self, seed):
        rng = np.random.default_rng(seed)
        adj_start, adj_stop = 1000, 1200
        n_frags = rng.integers(0, 200)

        # Fragments may start/stop outside the window (frag_array is
        # called with intersect_policy="any"), so sample from a padded
        # range and keep only those that actually overlap the window.
        starts = rng.integers(adj_start - 50, adj_stop + 50, size=n_frags)
        lengths = rng.integers(1, 300, size=n_frags)
        stops = starts + lengths
        strands = rng.integers(0, 2, size=n_frags).astype(bool)

        overlap_mask = (starts < adj_stop) & (stops > adj_start)
        starts, stops, strands = (
            starts[overlap_mask],
            stops[overlap_mask],
            strands[overlap_mask],
        )

        expected_depth, expected_ends = _brute_force_coverage_and_ends(
            starts, stops, strands, adj_start, adj_stop
        )
        depth, ends = _coverage_and_ends(
            starts, stops, strands, adj_start, adj_stop
        )

        np.testing.assert_array_equal(depth, expected_depth)
        np.testing.assert_array_equal(ends, expected_ends)

    def test_empty_fragments(self):
        starts = np.array([], dtype=np.int64)
        stops = np.array([], dtype=np.int64)
        strands = np.array([], dtype=bool)

        depth, ends = _coverage_and_ends(starts, stops, strands, 100, 150)

        assert np.all(depth == 0)
        assert np.all(ends == 0)
        assert depth.shape == (50,)
        assert ends.shape == (50,)

    def test_degenerate_interval(self):
        starts = np.array([100, 105], dtype=np.int64)
        stops = np.array([120, 130], dtype=np.int64)
        strands = np.array([True, False], dtype=bool)

        depth, ends = _coverage_and_ends(starts, stops, strands, 100, 100)

        assert depth.shape == (0,)
        assert ends.shape == (0,)

    def test_fragment_spans_entire_window(self):
        # A fragment starting before and ending after the window should
        # contribute to depth everywhere but not to ends anywhere.
        starts = np.array([50], dtype=np.int64)
        stops = np.array([200], dtype=np.int64)
        strands = np.array([True], dtype=bool)

        depth, ends = _coverage_and_ends(starts, stops, strands, 100, 150)

        assert np.all(depth == 1)
        assert np.all(ends == 0)

    def test_ends_exactly_at_window_boundaries(self):
        # Window is [100, 150). frag 0 (+): starts exactly at adj_start
        # (100) -> should register. frag 1 (-): stops exactly at adj_stop
        # (150), one past the last position (149) -> half-open, should
        # NOT register. frag 2 (-): stops at 149, the last position in
        # the window -> should register.
        starts = np.array([100, 90, 90], dtype=np.int64)
        stops = np.array([200, 150, 149], dtype=np.int64)
        strands = np.array([True, False, False], dtype=bool)

        depth, ends = _coverage_and_ends(starts, stops, strands, 100, 150)

        expected_depth, expected_ends = _brute_force_coverage_and_ends(
            starts, stops, strands, 100, 150
        )
        np.testing.assert_array_equal(depth, expected_depth)
        np.testing.assert_array_equal(ends, expected_ends)
        assert ends[0] == 1  # start of frag 0 at position 100
        assert ends[-1] == 1  # stop of frag 2 at position 149

