"""Cleavage Profiler

This module is used to find the cleavage profile as described in Zhou
et al 2022 (https://doi.org/10.1073/pnas.2209852119). Cleavage profile
describes the proportion of fragment ends at a site over the depth at
the site (cleavage proportion) calculated over a 5+/- nt window around a
CpG site.
"""

from __future__ import annotations
from typing import Union


def cleavage_profile(
    input_file: str,
    verbose: Union[int, bool]=False
):
    
    return NotImplementedError