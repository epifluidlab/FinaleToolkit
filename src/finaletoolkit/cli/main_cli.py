#!/usr/bin/env python3
"""
Entry point for the ``finaletoolkit`` command-line program.

The parser is assembled from the per-subcommand builders in
:mod:`finaletoolkit.cli.commands`.  Each subcommand records the implementing
``module``/``func`` via ``set_defaults``; :func:`main_cli` lazily imports that
function and calls it with the matching parsed arguments (filtered to the
function's signature), so importing the CLI stays cheap.
"""
from __future__ import annotations

import argparse
import importlib
from inspect import getfullargspec
from sys import stderr

import pysam

from finaletoolkit.io.reference import ReferenceWrapper
from finaletoolkit.utils.validation import validate_compatible_contigs

from .. import __version__
from .commands import build_subparsers


def main_cli_parser() -> argparse.ArgumentParser:
    """Build and return the top-level argument parser."""
    parser = argparse.ArgumentParser(
        description="FinaleToolkit is a package and standalone program to "
        "extract fragmentation features of cell-free DNA from paired-end "
        "sequencing data.",
        epilog="",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"FinaleToolkit {__version__}",
    )
    subparsers = parser.add_subparsers()
    build_subparsers(subparsers)
    return parser


def _validate_inputs(args: argparse.Namespace) -> None:
    """Validate CRAM reference requirements and contig compatibility."""
    input_file = getattr(args, "input_file", None)
    reference_file = getattr(
        args, "reference_file", getattr(args, "refseq_file", None)
    )

    if not input_file:
        return

    # CRAM input requires a reference.
    if str(input_file).lower().endswith(".cram") and not reference_file:
        stderr.write(
            "Error: CRAM files require a reference file (-r/--reference-file).\n"
        )
        exit(1)

    # When a reference is supplied for BAM/CRAM, check contig compatibility.
    if reference_file and (
        str(input_file).lower().endswith(".bam")
        or str(input_file).lower().endswith(".cram")
    ):
        try:
            try:
                with pysam.AlignmentFile(
                    input_file, "r", reference_filename=reference_file
                ) as sam:
                    input_contigs = dict(zip(sam.references, sam.lengths))
            except Exception as e:
                stderr.write(
                    f"Error opening alignment file '{input_file}': {e}\n"
                )
                exit(1)

            try:
                with ReferenceWrapper(reference_file) as ref:
                    ref_contigs = ref.chroms
            except Exception as e:
                stderr.write(
                    f"Error opening reference file '{reference_file}': {e}\n"
                )
                exit(1)

            validate_compatible_contigs(
                ref_contigs,
                input_contigs,
                validate_sizes=True,
                throw_on_error=True,
            )
        except ImportError:
            pass
        except (ValueError, RuntimeError) as e:
            stderr.write(f"Validation Error: {e}\n")
            exit(1)


def main_cli() -> None:
    """Parse arguments, validate inputs, and dispatch to the chosen subcommand."""
    parser = main_cli_parser()

    args = parser.parse_args()
    _validate_inputs(args)
    if hasattr(args, "func"):
        try:
            funcargs = vars(args)
            func_module = funcargs.pop("module")
            func_name = funcargs.pop("func")

            module = importlib.import_module(func_module)
            function = getattr(module, func_name)
            spec = getfullargspec(function)
            if spec.varkw is None:
                funcargs = {
                    key: value
                    for key, value in funcargs.items()
                    if key in spec.args or key in spec.kwonlyargs
                }

            function(**funcargs)
        except AttributeError as e:
            stderr.write(f"FinaleToolkit recieved AttributeError: {e}\n")
            stderr.write("Please see usage instructions below.\n")
            parser.print_help()
    else:
        parser.print_help()


if __name__ == "__main__":
    main_cli()
