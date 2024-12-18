"""
A package and standalone program to extract fragmentation features of
cell-free DNA from paired-end sequencing data.
"""

import lazy_loader as lazy

from .version import __version__


subpackages = ["cli", "frag", "genome", "utils"]

__getattr__, __dir__, _ = lazy.attach(__name__, subpackages)