#!/usr/bin/env python3
"""
CHANGELOG.md helpers used by the release-automation workflows.

  bump <version>
      Move the "## [Unreleased]" section's content under a new
      "## [<version>] - <date>" heading, leaving a fresh, empty
      "## [Unreleased]" section above it. Fails if [Unreleased] has no
      entries -- there's nothing to release.

  extract <version>
      Print the body of the "## [<version>]" section, for use as a
      GitHub Release's notes.

Used by ``.github/workflows/release-prepare.yml`` (bump, when preparing a
release PR) and ``.github/workflows/release-tag.yml`` (extract, when writing
the GitHub Release body after that PR merges).
"""
from __future__ import annotations

import argparse
import datetime
import pathlib
import re
import sys

CHANGELOG = pathlib.Path(__file__).resolve().parent.parent / "CHANGELOG.md"

# Matches a Keep a Changelog version heading, e.g. "## [Unreleased]" or
# "## [1.0.0] - 2026-06-26".
_HEADING_RE = re.compile(r"^## \[(?P<version>[^\]]+)\](?: - .+)?[ \t]*$", re.MULTILINE)


def _sections(text: str) -> list[dict]:
    """Return each "## [...]" heading's version and content span."""
    matches = list(_HEADING_RE.finditer(text))
    sections = []
    for i, m in enumerate(matches):
        content_end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        sections.append(
            {
                "version": m.group("version"),
                "heading_start": m.start(),
                "content_start": m.end(),
                "content_end": content_end,
            }
        )
    return sections


def bump(version: str) -> None:
    """Move [Unreleased]'s content under a new dated "## [version]" heading."""
    text = CHANGELOG.read_text()
    sections = _sections(text)
    if not sections or sections[0]["version"] != "Unreleased":
        sys.exit("CHANGELOG.md must start with an '## [Unreleased]' section.")

    unreleased = sections[0]
    body = text[unreleased["content_start"] : unreleased["content_end"]].strip("\n")
    if not body:
        sys.exit(
            "## [Unreleased] has no entries -- nothing to release. "
            "Add changelog entries before preparing a release."
        )

    if any(s["version"] == version for s in sections):
        sys.exit(f"CHANGELOG.md already has a '## [{version}]' section.")

    today = datetime.date.today().isoformat()
    replacement = f"## [Unreleased]\n\n## [{version}] - {today}\n\n{body}\n\n"
    updated = (
        text[: unreleased["heading_start"]]
        + replacement
        + text[unreleased["content_end"] :]
    )
    CHANGELOG.write_text(updated)
    print(f"Prepared CHANGELOG.md for {version} ({today}).")


def extract(version: str) -> None:
    """Print the body of the "## [version]" section."""
    text = CHANGELOG.read_text()
    for section in _sections(text):
        if section["version"] == version:
            body = text[section["content_start"] : section["content_end"]]
            print(body.strip("\n"))
            return
    sys.exit(f"No '## [{version}]' section found in CHANGELOG.md.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    bump_parser = subparsers.add_parser(
        "bump", help="Move [Unreleased] content under a new dated version heading."
    )
    bump_parser.add_argument(
        "version", help='Version being released, e.g. "1.1.0" (no leading "v").'
    )

    extract_parser = subparsers.add_parser(
        "extract", help="Print a version's changelog section (for release notes)."
    )
    extract_parser.add_argument("version")

    args = parser.parse_args()
    if args.command == "bump":
        bump(args.version)
    elif args.command == "extract":
        extract(args.version)


if __name__ == "__main__":
    main()
