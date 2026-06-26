"""
Command-line interface for FinaleToolkit.

``main_cli`` (the Click command group) is resolved lazily so that ``import
finaletoolkit.cli`` does not eagerly import the ``main_cli`` submodule (which
keeps ``python -m finaletoolkit.cli.main_cli`` free of double-import warnings).
"""
from __future__ import annotations

__all__ = ["main_cli"]


def __getattr__(name: str):
    if name in __all__:
        import importlib

        # Import by full dotted path (not ``from . import main_cli``) so the
        # submodule lookup cannot re-enter this ``__getattr__``. The submodule
        # shares the name ``main_cli`` with the attribute we expose.
        module = importlib.import_module(f"{__name__}.main_cli")
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
