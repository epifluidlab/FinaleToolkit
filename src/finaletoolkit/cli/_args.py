"""
Reusable argument groups for the FinaleToolkit CLI.

This is a redesigned, consistent flag scheme (a clean break from the original
mixed naming). Conventions applied uniformly across every subcommand:

* ``-o/--output``        output path (was ``--output-file``)
* ``-r/--reference``     optional reference genome (was ``--reference-file``)
* ``-q/--min-mapq``      minimum mapping quality (was ``-q/--quality-threshold``)
* ``--min-length`` / ``--max-length``   fragment-length bounds (no more
  ``-min``/``-max``/``-lo``/``-hi``/``--fraction-*``)
* ``-t/--threads``       worker processes (was ``-w/--workers``)
* ``-v/--verbose``       counting verbosity (was a mix of store_true/count)
* ``-k/--kmer-length``   k-mer length (was bare ``-k``)
* ``--strand {both,forward,reverse}``  strand selection for motif commands
  (replaces the confusing ``--single-strand``/``--negative-strand`` toggle)

Boolean on/off options use :class:`argparse.BooleanOptionalAction` so the
``--x``/``--no-x`` pair shows a sensible default, instead of a lone negated flag
whose ``store_false`` ``dest`` rendered a misleading ``Default: True``.

Each flag keeps the ``dest`` that matches its Python function parameter, so the
Python API is unchanged — only the command-line spelling changed.
"""
from __future__ import annotations

import argparse

__all__ = [
    "add_input_file",
    "add_reference_option",
    "add_output",
    "add_min_mapq",
    "add_min_length",
    "add_max_length",
    "add_intersect_policy",
    "add_threads",
    "add_verbose",
    "add_kmer",
    "add_strand",
    "add_bool_flag",
]

_INPUT_HELP = "Path to a BAM/CRAM/Fragment file containing fragment data."
_REFERENCE_HELP = (
    "Path to a FASTA (.fa, .fasta, .fna) reference genome file. "
    "Required for CRAM input."
)


def add_input_file(parser: argparse.ArgumentParser, help: str = _INPUT_HELP) -> None:
    """Add the positional ``input_file`` argument."""
    parser.add_argument("input_file", metavar="INPUT", help=help)


def add_reference_option(parser: argparse.ArgumentParser) -> None:
    """Add the optional ``-r/--reference`` argument (dest ``reference_file``)."""
    parser.add_argument(
        "-r",
        "--reference",
        dest="reference_file",
        metavar="FASTA",
        required=False,
        help=_REFERENCE_HELP,
    )


def add_output(parser: argparse.ArgumentParser, help: str, default: str = "-") -> None:
    """Add the ``-o/--output`` argument (dest ``output_file``)."""
    parser.add_argument(
        "-o",
        "--output",
        dest="output_file",
        metavar="PATH",
        default=default,
        help=help,
    )


def add_min_mapq(parser: argparse.ArgumentParser, default: int = 30) -> None:
    """Add ``-q/--min-mapq`` (dest ``quality_threshold``)."""
    parser.add_argument(
        "-q",
        "--min-mapq",
        dest="quality_threshold",
        metavar="MAPQ",
        default=default,
        type=int,
        help="Minimum mapping quality (MAPQ) for a fragment to be included.",
    )


def add_min_length(parser: argparse.ArgumentParser, default, help: str) -> None:
    """Add ``--min-length`` (dest ``min_length``)."""
    parser.add_argument(
        "--min-length", dest="min_length", metavar="BP", default=default, type=int,
        help=help,
    )


def add_max_length(parser: argparse.ArgumentParser, default, help: str) -> None:
    """Add ``--max-length`` (dest ``max_length``)."""
    parser.add_argument(
        "--max-length", dest="max_length", metavar="BP", default=default, type=int,
        help=help,
    )


def add_intersect_policy(parser: argparse.ArgumentParser) -> None:
    """Add ``-p/--intersect-policy`` (choices midpoint/any)."""
    parser.add_argument(
        "-p",
        "--intersect-policy",
        choices=["midpoint", "any"],
        default="midpoint",
        type=str,
        help="How a fragment is counted as inside an interval: 'midpoint' (the "
        "fragment midpoint lies in the interval) or 'any' (any overlap). "
        "See the User Guide for details.",
    )


def add_threads(parser: argparse.ArgumentParser, default: int = 1) -> None:
    """Add ``-t/--threads`` (dest ``workers``)."""
    parser.add_argument(
        "-t",
        "--threads",
        dest="workers",
        metavar="N",
        default=default,
        type=int,
        help="Number of worker processes to use.",
    )


def add_verbose(parser: argparse.ArgumentParser) -> None:
    """Add the counting ``-v/--verbose`` argument."""
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity (repeatable, e.g. -vv) for detailed processing "
        "information.",
    )


def add_kmer(parser: argparse.ArgumentParser, default: int) -> None:
    """Add ``-k/--kmer-length`` (dest ``k``)."""
    parser.add_argument(
        "-k",
        "--kmer-length",
        dest="k",
        metavar="K",
        default=default,
        type=int,
        help="Length of the k-mer motif.",
    )


class _StrandAction(argparse.Action):
    """Translate ``--strand {both,forward,reverse}`` into the two function args.

    The motif functions take ``both_strands`` and ``negative_strand`` booleans;
    exposing those directly on the CLI as a ``store_false`` toggle produced a
    confusing ``Default: True``.  This action maps a single, explicit choice
    onto both destinations instead.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, "both_strands", values == "both")
        setattr(namespace, "negative_strand", values == "reverse")


def add_strand(parser: argparse.ArgumentParser, feature: str) -> None:
    """Add the ``--strand {both,forward,reverse}`` selector.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The subparser to extend.
    feature : str
        Wording for help text, e.g. ``"end"`` or ``"breakpoint"``.
    """
    parser.add_argument(
        "--strand",
        choices=["both", "forward", "reverse"],
        default="both",
        action=_StrandAction,
        help=f"Which fragment strand(s) to use for {feature} motifs: 'both' "
        "(default), 'forward' (5' ends of the positive strand only), or "
        "'reverse' (5' ends of the negative strand only).",
    )
    # Defaults for the destinations the action writes to (i.e. 'both').
    parser.set_defaults(both_strands=True, negative_strand=False)


def add_bool_flag(
    parser: argparse.ArgumentParser,
    name: str,
    dest: str,
    default: bool,
    help: str,
) -> None:
    """Add a ``--name`` / ``--no-name`` boolean flag.

    Uses :class:`argparse.BooleanOptionalAction` so the default is reported
    against the positive form (``--name``) and reads naturally, rather than the
    confusing ``Default: True`` produced by a lone ``store_false`` flag.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The subparser to extend.
    name : str
        The positive flag name without leading dashes (e.g. ``"savgol"``).
    dest : str
        The destination (must match the function parameter name).
    default : bool
        Default value.
    help : str
        Help text describing the positive behavior.
    """
    parser.add_argument(
        f"--{name}",
        dest=dest,
        default=default,
        action=argparse.BooleanOptionalAction,
        help=help,
    )
