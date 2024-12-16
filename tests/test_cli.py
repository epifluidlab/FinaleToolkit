"""
Tests for main_cli and entry points.
"""
import os
from inspect import getfullargspec

import pytest

from finaletoolkit.cli.main_cli import main_cli_parser

IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"

parser = main_cli_parser()
subcommands: dict = parser._subparsers._actions[2]._name_parser_map

class TestCLIArgs:
    """
    Test if provided commandline flags match args in associated function.
    """
    @pytest.mark.skipif(
        IN_GITHUB_ACTIONS,
        reason="Test doesn't always work in Github Actions.")
    @pytest.mark.parametrize("name,subparser", subcommands.items())
    def test_cli_args(self, name, subparser):
        # find args and associated func for each subparser
        # get args for CLI
        cli_args = [action.dest for action in subparser._actions[1:]]

        # get args for func
        func = subparser._defaults['func']
        func_args = getfullargspec(func).args

        # check cli args are subset of func args
        for arg in cli_args:
            assert arg in func_args, f"CLI arg {arg} of {name} not in {func}"
            
        # check func args are subset of cli args
        for arg in func_args:
            assert arg in cli_args or arg == 'fraction_low' or arg == 'fraction_high', f"API arg {arg} of {func} not in {name}"


class TestCLIEntryPoint:
    """
    Test all CLI subcommands related to end_motifs and MDS, genomewide
    and intervals.
    """

    def test_coverage(self):
        exit_status = os.system('finaletoolkit coverage --help')
        assert exit_status == 0

    def test_frag_length_bins(self):
        exit_status = os.system('finaletoolkit frag-length-bins --help')
        assert exit_status == 0

    def test_frag_length_intervals(self):
        exit_status = os.system('finaletoolkit frag-length-intervals --help')
        assert exit_status == 0

    def test_cleavage_profile(self):
        exit_status = os.system('finaletoolkit cleavage-profile --help')
        assert exit_status == 0

    def test_wps(self):
        exit_status = os.system('finaletoolkit wps --help')
        assert exit_status == 0

    def test_agg_wps(self):
        exit_status = os.system('finaletoolkit adjust-wps --help')
        assert exit_status == 0

    def test_delfi(self):
        exit_status = os.system('finaletoolkit delfi --help')
        assert exit_status == 0

    def test_delfi_gc_correct(self):
        exit_status = os.system('finaletoolkit delfi-gc-correct --help')
        assert exit_status == 0

    def test_end_motif(self):
        exit_status = os.system('finaletoolkit end-motifs --help')
        assert exit_status == 0

    def test_mds(self):
        exit_status = os.system('finaletoolkit mds --help')
        assert exit_status == 0

    def test_interval_end_motif(self):
        exit_status = os.system('finaletoolkit interval-end-motifs --help')
        assert exit_status == 0

    def test_interval_mds(self):
        exit_status = os.system('finaletoolkit interval-mds --help')
        assert exit_status == 0

    def test_filter_bam(self):
        exit_status = os.system('finaletoolkit filter-bam --help')
        assert exit_status == 0

    def test_agg_bw(self):
        exit_status = os.system('finaletoolkit agg-bw --help')
        assert exit_status == 0

    def test_gap_bed(self):
        exit_status = os.system('finaletoolkit gap-bed --help')
        assert exit_status == 0
