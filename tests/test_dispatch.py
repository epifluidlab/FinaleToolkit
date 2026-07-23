"""
Tests for finaletoolkit.cli._dispatch

Covers the pure, CRAM-independent logic: strand-flag translation, the
input-file-absent/BAM-without-reference fast paths through
``_validate_inputs``, and ``run``'s param-filtering/dispatch mechanics. The
CRAM+reference contig-validation branch (which needs real alignment/reference
files) is already exercised by the filter-file tests in test_cli.py.
"""

import pytest

from finaletoolkit.cli._dispatch import _translate_strand, _validate_inputs, run


class TestTranslateStrand:
    def test_both(self):
        params = {"strand": "both"}
        _translate_strand(params)
        assert params == {"both_strands": True, "negative_strand": False}

    def test_forward(self):
        params = {"strand": "forward"}
        _translate_strand(params)
        assert params == {"both_strands": False, "negative_strand": False}

    def test_reverse(self):
        params = {"strand": "reverse"}
        _translate_strand(params)
        assert params == {"both_strands": False, "negative_strand": True}

    def test_no_strand_key_is_a_no_op(self):
        params = {"other": 1}
        _translate_strand(params)
        assert params == {"other": 1}


class TestValidateInputs:
    def test_no_input_file_is_a_no_op(self):
        _validate_inputs({})
        _validate_inputs({"input_file": None})
        _validate_inputs({"input_file": ""})

    def test_cram_without_reference_exits(self):
        with pytest.raises(SystemExit) as exc_info:
            _validate_inputs({"input_file": "sample.cram"})
        assert exc_info.value.code == 1

    def test_bam_without_reference_is_a_no_op(self):
        # A reference is optional for BAM; only CRAM requires one.
        _validate_inputs({"input_file": "sample.bam"})

    def test_non_alignment_input_is_a_no_op(self):
        # Neither .bam nor .cram: the reference/contig-compatibility branch
        # doesn't apply regardless of whether a reference was given.
        _validate_inputs({"input_file": "sample.frag.gz", "reference_file": "ref.fa"})


def _stub_add(a, b=10):
    return a + b


def _stub_kwargs(**kwargs):
    return kwargs


class TestRun:
    def test_filters_params_to_function_signature(self):
        result = run(
            __name__, "_stub_add", {"a": 1, "b": 2, "unrelated_cli_only_key": "x"}
        )
        assert result == 3

    def test_uses_function_defaults_for_missing_params(self):
        result = run(__name__, "_stub_add", {"a": 5})
        assert result == 15

    def test_varkw_function_receives_all_params_unfiltered(self):
        result = run(__name__, "_stub_kwargs", {"a": 1, "anything": "goes"})
        assert result == {"a": 1, "anything": "goes"}

    def test_strand_translated_before_dispatch(self):
        result = run(
            __name__, "_stub_kwargs", {"strand": "reverse", "a": 1}
        )
        assert result == {"a": 1, "both_strands": False, "negative_strand": True}
