"""
Command-line interface for FinaleToolkit.

``main_cli`` and ``main_cli_parser`` are imported lazily so that
``python -m finaletoolkit.cli.main_cli`` does not trigger a double-import
RuntimeWarning.
"""
from __future__ import annotations

__all__ = ["main_cli", "main_cli_parser"]


def __getattr__(name: str):
    if name in __all__:
        from . import main_cli as _module

        return getattr(_module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
