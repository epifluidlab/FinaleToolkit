"""
Tests for main_cli and entry points.
"""
from __future__ import annotations
import os
from pathlib import Path
from inspect import getfullargspec
import importlib
import subprocess
import sys

import pytest

from finaletoolkit.cli.commands import COMMAND_TARGETS, COMMANDS

# name -> Click command object, for introspecting declared parameters.
_COMMANDS_BY_NAME = {command.name: command for command in COMMANDS}

class TestCLIArgs:
    """
    Test if provided commandline flags match args in associated function.
    """
    @pytest.mark.parametrize("name,target", COMMAND_TARGETS.items())
    def test_lazy_import(self, name: str, target: tuple[str, str]):
        # the (module, func) the subcommand dispatches to
        module_path, func = target

        # try to see module spec
        module = importlib.import_module(module_path)
        function = getattr(module, func)
        assert callable(function), f'The {func} is not a callable in {module}.'

    # Function arguments intentionally not exposed as CLI flags:
    # fraction_low/fraction_high are deprecated aliases of the length filters,
    # and gc_correct is the deprecated alias of delfi's --no-gc-correct.
    _UNEXPOSED_ARGS = {"fraction_low", "fraction_high", "gc_correct"}

    @pytest.mark.parametrize("name,target", COMMAND_TARGETS.items())
    def test_cli_args(self, name: str, target: tuple[str, str]):
        """Every CLI flag maps to a real function argument, and every function
        argument is exposed as a flag (modulo the deprecated/internal ones).

        Guards against a misrouted Click ``dest`` silently passing the wrong
        value, or a feature parameter being dropped from the CLI.  Commands
        whose target is a ``_cli_*`` wrapper (mds, regional-mds, gap-bed) are
        CLI-specific shims, not 1:1 with a feature function, so they are only
        checked for importability (see ``test_lazy_import``).
        """
        module_path, func_name = target
        func = getattr(importlib.import_module(module_path), func_name)

        if func_name.startswith("_cli_"):
            assert callable(func), f"{name} target {func_name} is not callable"
            return

        # Click parameter names equal the function argument names; --strand
        # expands to both_strands/negative_strand at dispatch time.
        command = _COMMANDS_BY_NAME[name]
        cli_args = [param.name for param in command.params]
        if "strand" in cli_args:
            cli_args.remove("strand")
            cli_args += ["both_strands", "negative_strand"]

        spec = getfullargspec(func)
        func_args = set(spec.args) | set(spec.kwonlyargs)

        # forward: no CLI flag may reference a non-existent function argument
        for arg in cli_args:
            assert arg in func_args, f"CLI arg {arg!r} of {name} not in {func}"

        # reverse: no function argument may be silently dropped from the CLI
        for arg in func_args:
            assert arg in cli_args or arg in self._UNEXPOSED_ARGS, (
                f"API arg {arg!r} of {func} not exposed by {name}")


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

    def test_end_motif(self):
        exit_status = os.system('finaletoolkit end-motifs --help')
        assert exit_status == 0

    def test_mds(self):
        exit_status = os.system('finaletoolkit mds --help')
        assert exit_status == 0

    def test_interval_end_motif(self):
        exit_status = os.system('finaletoolkit interval-end-motifs --help')
        assert exit_status == 0

    def test_regional_mds(self):
        exit_status = os.system('finaletoolkit regional-mds --help')
        assert exit_status == 0

    def test_filter_file(self):
        exit_status = os.system('finaletoolkit filter-file --help')
        assert exit_status == 0

    def test_agg_bw(self):
        exit_status = os.system('finaletoolkit agg-bw --help')
        assert exit_status == 0

    def test_gap_bed(self):
        exit_status = os.system('finaletoolkit gap-bed --help')
        assert exit_status == 0

    def test_coverage_smoke(self, request):
        root = Path(request.config.rootpath)
        input_file = root / 'tests' / 'data' / '12.3444.b37.frag.gz'
        interval_file = root / 'tests' / 'data' / 'intervals.bed'
        result = subprocess.run(
            [
                sys.executable,
                '-m',
                'finaletoolkit.cli.main_cli',
                'coverage',
                str(input_file),
                str(interval_file),
                '--normalize',
                '-o',
                '-',
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stderr
        assert result.stdout.splitlines() == [
            '12\t34443118\t34443538\t.\t0.25',
            '12\t34444968\t34446115\t.\t0.4375',
        ]
