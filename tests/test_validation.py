"""
Tests for finaletoolkit.utils.validation
"""

import pytest

from finaletoolkit.utils.validation import (
    validate_compatible_contigs,
    valid_interval,
)


class TestValidateCompatibleContigs:
    def test_exact_match_lists(self):
        assert validate_compatible_contigs(["1", "2", "3"], ["1", "2", "3"])

    def test_subset_allowed_by_default(self):
        assert validate_compatible_contigs(["1", "2", "3"], ["1", "2"])

    def test_extra_input_contig_raises_by_default(self):
        with pytest.raises(ValueError, match="not found in reference"):
            validate_compatible_contigs(["1", "2"], ["1", "2", "4"])

    def test_extra_input_contig_returns_false_when_not_throwing(self):
        assert not validate_compatible_contigs(
            ["1", "2"], ["1", "2", "4"], throw_on_error=False
        )

    def test_subset_disallowed_raises(self):
        with pytest.raises(ValueError, match="not found in input"):
            validate_compatible_contigs(
                ["1", "2", "3"], ["1", "2"], allow_subset=False
            )

    def test_subset_disallowed_returns_false_when_not_throwing(self):
        assert not validate_compatible_contigs(
            ["1", "2", "3"], ["1", "2"], allow_subset=False, throw_on_error=False
        )

    def test_subset_disallowed_but_sets_match_passes(self):
        assert validate_compatible_contigs(
            ["1", "2"], ["1", "2"], allow_subset=False
        )

    def test_validate_sizes_matching(self):
        assert validate_compatible_contigs(
            {"1": 100, "2": 200}, {"1": 100, "2": 200}, validate_sizes=True
        )

    def test_validate_sizes_mismatch_raises(self):
        with pytest.raises(RuntimeError, match="length mismatch"):
            validate_compatible_contigs(
                {"1": 100, "2": 200}, {"1": 100, "2": 999}, validate_sizes=True
            )

    def test_validate_sizes_mismatch_returns_false_when_not_throwing(self):
        assert not validate_compatible_contigs(
            {"1": 100, "2": 200},
            {"1": 100, "2": 999},
            validate_sizes=True,
            throw_on_error=False,
        )

    def test_validate_sizes_requires_dicts_raises_typeerror(self):
        with pytest.raises(TypeError, match="requires both"):
            validate_compatible_contigs(
                ["1", "2"], ["1", "2"], validate_sizes=True
            )

    def test_validate_sizes_requires_dicts_returns_false_when_not_throwing(self):
        assert not validate_compatible_contigs(
            ["1", "2"], ["1", "2"], validate_sizes=True, throw_on_error=False
        )

    def test_validate_sizes_only_checks_input_contigs(self):
        # Reference may have extra contigs (allowed by default subset rule);
        # their sizes shouldn't be checked since they aren't in input_names.
        assert validate_compatible_contigs(
            {"1": 100, "2": 999999},
            {"1": 100},
            validate_sizes=True,
        )


class TestValidInterval:
    def test_contig_missing_returns_false_by_default(self):
        assert not valid_interval(["1", "2"], "3")

    def test_contig_missing_raises_when_requested(self):
        with pytest.raises(ValueError, match="not found in reference"):
            valid_interval(["1", "2"], "3", throw_on_error=True)

    def test_valid_bounds_with_lengths(self):
        assert valid_interval({"1": 1000}, "1", start=0, stop=1000)

    def test_start_negative_with_lengths_returns_false(self):
        assert not valid_interval({"1": 1000}, "1", start=-1)

    def test_start_negative_with_lengths_raises(self):
        with pytest.raises(IndexError, match="out of bounds"):
            valid_interval({"1": 1000}, "1", start=-1, throw_on_error=True)

    def test_start_beyond_length_returns_false(self):
        assert not valid_interval({"1": 1000}, "1", start=1000)

    def test_stop_negative_returns_false(self):
        assert not valid_interval({"1": 1000}, "1", stop=-1)

    def test_stop_beyond_length_raises(self):
        with pytest.raises(IndexError, match="out of bounds"):
            valid_interval({"1": 1000}, "1", stop=1001, throw_on_error=True)

    def test_start_not_less_than_stop_returns_false(self):
        assert not valid_interval({"1": 1000}, "1", start=500, stop=500)

    def test_start_not_less_than_stop_raises(self):
        with pytest.raises(ValueError, match="must be less than stop"):
            valid_interval(
                {"1": 1000}, "1", start=500, stop=100, throw_on_error=True
            )

    def test_only_start_provided_skips_stop_checks(self):
        assert valid_interval({"1": 1000}, "1", start=999)

    def test_only_stop_provided_skips_start_checks(self):
        assert valid_interval({"1": 1000}, "1", stop=1000)

    def test_list_input_contig_present_no_length_checks(self):
        # Without lengths, only a negative start can be rejected.
        assert valid_interval(["1", "2"], "1", start=10**9)

    def test_list_input_negative_start_returns_false(self):
        assert not valid_interval(["1", "2"], "1", start=-1)

    def test_list_input_negative_start_raises(self):
        with pytest.raises(IndexError, match="cannot be negative"):
            valid_interval(["1", "2"], "1", start=-1, throw_on_error=True)

    def test_list_input_no_start_or_stop_is_valid(self):
        assert valid_interval(["1", "2"], "1")
