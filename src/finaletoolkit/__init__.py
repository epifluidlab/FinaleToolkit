"""
A package and standalone program to extract fragmentation features of
cell-free DNA from paired-end sequencing data.
"""

from .version import __version__


__all__ = ["cli", "frag", "genome", "utils"]

# Delay imports until the submodule is actually accessed.
def __getattr__(name):
    if name in __all__:
        import importlib
        return importlib.import_module(f".{name}", __name__)
    raise AttributeError(f"Module {__name__} has no attribute {name}")
