"""
Single-source module for the package version number.

The version string is derived from git tags/commit distance by
setuptools-scm at build/install time, written to the generated
``finaletoolkit._version`` module. That file doesn't exist in a source
checkout that hasn't been built or installed yet (e.g. running from a
raw ``git clone`` without ``pip install``), hence the fallback below.
"""

try:
    from ._version import __version__
except ImportError:
    __version__ = "0.0.0.dev0+unknown"
