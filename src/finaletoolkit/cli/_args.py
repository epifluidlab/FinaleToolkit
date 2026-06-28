"""
Reusable Click options/arguments for the FinaleToolkit CLI.

Each helper returns a Click decorator so the per-command definitions in
:mod:`finaletoolkit.cli.commands` stay short and the flag scheme stays uniform
across every subcommand:

* ``-o/--output``        output path        (param ``output_file``)
* ``-r/--reference``     optional reference (param ``reference_file``)
* ``-q/--min-mapq``      minimum MAPQ       (param ``quality_threshold``)
* ``--min-length`` / ``--max-length``       fragment-length bounds
* ``-t/--threads``       worker processes   (param ``workers``)
* ``-v/--verbose``       counting verbosity (param ``verbose``)
* ``-k/--kmer-length``   k-mer length       (param ``k``)
* ``--strand {both,forward,reverse}``       strand selection for motif commands

The first string after the flag spellings is the *parameter name*, which is
deliberately kept equal to the implementing function's argument.  Only the
command-line spelling is new; the Python API is unchanged.
"""
from __future__ import annotations

import rich_click as click

__all__ = [
    "input_argument",
    "reference_option",
    "output_option",
    "min_mapq_option",
    "min_length_option",
    "max_length_option",
    "intersect_policy_option",
    "threads_option",
    "verbose_option",
    "kmer_option",
    "strand_option",
    "bool_flag",
]

_INPUT_HELP = "Path to a BAM/CRAM/Fragment file containing fragment data."
_REFERENCE_HELP = (
    "FASTA (.fa, .fasta, .fna) reference genome. Required for CRAM input."
)


def input_argument(metavar: str = "INPUT"):
    """Positional ``input_file`` (param ``input_file``)."""
    return click.argument("input_file", metavar=metavar)


def reference_option():
    """Optional ``-r/--reference`` (param ``reference_file``)."""
    return click.option(
        "-r",
        "--reference",
        "reference_file",
        metavar="FASTA",
        default=None,
        help=_REFERENCE_HELP,
    )


def output_option(help: str, default: str = "-"):
    """``-o/--output`` (param ``output_file``).

    The shared ``-`` convention (write to stdout) is appended to every output
    help string so it is documented consistently across commands.
    """
    return click.option(
        "-o",
        "--output",
        "output_file",
        metavar="PATH",
        default=default,
        show_default=True,
        help=f"{help} Pass '-' to write to standard output (stdout).",
    )


def min_mapq_option(default: int = 30):
    """``-q/--min-mapq`` (param ``quality_threshold``)."""
    return click.option(
        "-q",
        "--min-mapq",
        "quality_threshold",
        metavar="MAPQ",
        default=default,
        show_default=True,
        type=int,
        help="Minimum mapping quality (MAPQ) for a fragment to be included.",
    )


def min_length_option(default, help: str):
    """``--min-length`` (param ``min_length``)."""
    return click.option(
        "--min-length",
        "min_length",
        metavar="BP",
        default=default,
        show_default=default is not None,
        type=int,
        help=help,
    )


def max_length_option(default, help: str):
    """``--max-length`` (param ``max_length``)."""
    return click.option(
        "--max-length",
        "max_length",
        metavar="BP",
        default=default,
        show_default=default is not None,
        type=int,
        help=help,
    )


def intersect_policy_option():
    """``-p/--intersect-policy`` (param ``intersect_policy``)."""
    return click.option(
        "-p",
        "--intersect-policy",
        "intersect_policy",
        type=click.Choice(["midpoint", "any"]),
        default="midpoint",
        show_default=True,
        help="How a fragment counts as inside an interval: 'midpoint' (the "
        "fragment midpoint lies in the interval) or 'any' (any overlap).",
    )


def threads_option(default: int = 1):
    """``-t/--threads`` (param ``workers``)."""
    return click.option(
        "-t",
        "--threads",
        "workers",
        metavar="N",
        default=default,
        show_default=True,
        type=int,
        help="Number of worker processes to use.",
    )


def verbose_option():
    """``-v/--verbose`` counting flag (param ``verbose``)."""
    return click.option(
        "-v",
        "--verbose",
        "verbose",
        count=True,
        default=0,
        metavar="",
        help="Increase verbosity (repeatable, e.g. -vv) for detailed output.",
    )


def kmer_option(default: int):
    """``-k/--kmer-length`` (param ``k``)."""
    return click.option(
        "-k",
        "--kmer-length",
        "k",
        metavar="K",
        default=default,
        show_default=True,
        type=int,
        help="Length of the k-mer motif.",
    )


def strand_option(feature: str):
    """``--strand {both,forward,reverse}`` selector for motif commands.

    Expanded into ``both_strands``/``negative_strand`` at dispatch time (see
    :func:`finaletoolkit.cli._dispatch._translate_strand`).
    """
    return click.option(
        "--strand",
        "strand",
        type=click.Choice(["both", "forward", "reverse"]),
        default="both",
        show_default=True,
        help=f"Fragment strand(s) for {feature} motifs: 'both', 'forward' (5' "
        "ends of the positive strand only), or 'reverse' (negative strand only).",
    )


def bool_flag(name: str, param: str, default: bool, help: str):
    """A ``--name/--no-name`` boolean flag (param ``param``).

    The default is stated in the help text ("On by default" / "Off by default")
    rather than shown as a bare ``True``/``False``, which is clearer in both the
    terminal and the generated docs.
    """
    state = "On" if default else "Off"
    return click.option(
        f"--{name}/--no-{name}",
        param,
        default=default,
        show_default=False,
        help=f"{help} {state} by default.",
    )
