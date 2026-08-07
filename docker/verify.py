"""Post-install smoke checks, run as a build step in the Dockerfile.

These catch packaging failures that otherwise surface only when a user runs
the subcommand that happens to need the missing piece. Kept as a file rather
than an inline `python -c` so it stays readable and can be run by hand:

    docker run --rm --entrypoint python <image> - < docker/verify.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import finaletoolkit

# Declared in [tool.setuptools.package-data]. If include-package-data or the
# package-data globs regress, the wheel installs fine and these vanish.
#
# Checked as plain paths under the installed package rather than via
# importlib.resources, which would import finaletoolkit.frag / .genome and
# pull in numpy, numba and friends -- turning a data-file check into a
# general import check. The CLI import below covers that separately.
REQUIRED_DATA = [
    "frag/data/end_motif_f_profiles.tsv",
    "genome/data/hg19.gap.txt.gz",
    "genome/data/hg38.gap.txt.gz",
]

failures: list[str] = []

version = finaletoolkit.__version__
print(f"finaletoolkit {version}")
if version.startswith("0.0.0.dev0"):
    failures.append(
        "version is the setuptools-scm fallback -- pass --build-arg "
        "SETUPTOOLS_SCM_PRETEND_VERSION=<version> (docker/build.sh does "
        "this for you)"
    )

package_root = Path(finaletoolkit.__file__).parent
for relative in REQUIRED_DATA:
    if not (package_root / relative).is_file():
        failures.append(f"missing package data: {relative}")

# The console script is what ENTRYPOINT invokes; importing the module it
# points at confirms the entry point resolves and that the compiled
# dependencies (numpy, pysam, numba, pyBigWig, py2bit) all load.
try:
    from finaletoolkit.cli.main_cli import main_cli  # noqa: F401
except Exception as exc:  # pragma: no cover - build-time guard
    failures.append(f"CLI entry point failed to import: {exc!r}")

if failures:
    print("\nbuild aborted:", file=sys.stderr)
    for failure in failures:
        print(f"  - {failure}", file=sys.stderr)
    sys.exit(1)

print(f"verified: {len(REQUIRED_DATA)} data files, CLI entry point importable")
