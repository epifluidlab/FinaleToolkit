#!/usr/bin/env python3
"""
Entry point for the ``finaletoolkit`` command-line program.

The CLI is built with `Click <https://click.palletsprojects.com>`_ and rendered
with `rich-click <https://github.com/ewels/rich-click>`_ for a polished help
display.  The subcommands live in :mod:`finaletoolkit.cli.commands`; each one
dispatches to its implementing function lazily (see
:mod:`finaletoolkit.cli._dispatch`), so ``--help`` stays fast.
"""
from __future__ import annotations

import rich_click as click

from .. import __version__
from .commands import register_commands

# --------------------------------------------------------------------------- #
# rich-click appearance: a restrained, consistent theme with accent colors for
# options/switches/metavars, soft panel borders, grouped command panels.
# --------------------------------------------------------------------------- #
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = False
click.rich_click.SHOW_METAVARS_COLUMN = True
click.rich_click.APPEND_METAVARS_HELP = False
click.rich_click.MAX_WIDTH = 100

# Monochrome help: no colors, so it stays legible on any terminal background
# (dark, light, or low-color). Structure comes from bold text and the panel
# boxes, not color. Every style that rich-click defaults to a color (cyan
# options/commands, yellow metavars/usage, red required markers) is overridden
# with a plain bold / dim / default style.
click.rich_click.STYLE_OPTION = "bold"
click.rich_click.STYLE_ARGUMENT = "bold"
click.rich_click.STYLE_COMMAND = "bold"
click.rich_click.STYLE_SWITCH = "bold"
click.rich_click.STYLE_METAVAR = "dim"
click.rich_click.STYLE_METAVAR_APPEND = "dim"
click.rich_click.STYLE_METAVAR_SEPARATOR = "dim"
click.rich_click.STYLE_USAGE = "bold"
click.rich_click.STYLE_USAGE_COMMAND = "bold"
click.rich_click.STYLE_HELPTEXT_FIRST_LINE = "bold"
click.rich_click.STYLE_HELPTEXT = ""
click.rich_click.STYLE_OPTION_DEFAULT = "dim"
click.rich_click.STYLE_OPTION_ENVVAR = "dim"
click.rich_click.STYLE_REQUIRED_SHORT = "bold"
click.rich_click.STYLE_REQUIRED_LONG = "dim"
click.rich_click.STYLE_OPTIONS_PANEL_BORDER = "dim"
click.rich_click.STYLE_COMMANDS_PANEL_BORDER = "dim"
click.rich_click.STYLE_OPTIONS_TABLE_BOX = "ROUNDED"
click.rich_click.STYLE_COMMANDS_TABLE_BOX = "SIMPLE"
click.rich_click.OPTIONS_PANEL_TITLE = "Options"
click.rich_click.ARGUMENTS_PANEL_TITLE = "Arguments"

# Organize the top-level command list into themed panels.
click.rich_click.COMMAND_GROUPS = {
    "finaletoolkit": [
        {
            "name": "Coverage & Fragment Length",
            "commands": [
                "coverage",
                "frag-length-bins",
                "frag-length-intervals",
            ],
        },
        {
            "name": "Protection & Cleavage",
            "commands": ["wps", "adjust-wps", "cleavage-profile"],
        },
        {
            "name": "DELFI",
            "commands": ["delfi"],
        },
        {
            "name": "Motifs & MDS",
            "commands": [
                "end-motifs",
                "interval-end-motifs",
                "breakpoint-motifs",
                "interval-breakpoint-motifs",
                "mds",
                "regional-mds",
            ],
        },
        {
            "name": "Utilities",
            "commands": ["filter-file", "agg-bw", "gap-bed"],
        },
    ]
}


@click.group(
    name="finaletoolkit",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.version_option(
    __version__,
    "-v",
    "--version",
    prog_name="FinaleToolkit",
    message="%(prog)s %(version)s",
)
def main_cli() -> None:
    """A package and standalone program to extract fragmentation features of
    cell-free DNA from paired-end sequencing data. Run a subcommand with
    --help to see its options and an example invocation. For file outputs,
    pass "-o -" to write to standard output instead of a file.
    """


register_commands(main_cli)


if __name__ == "__main__":
    main_cli()
