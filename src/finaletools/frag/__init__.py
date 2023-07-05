"""
Author: James Li
Created: 6/2/23
PI: Yaping Liu

Description:
Python script to calculate fragment features given a BAM file.

"""
# TODO: typing annotations for all functions

from __future__ import annotations
import argparse

from finaletools.frag.frag_length import frag_length
from finaletools.frag.coverage import coverage
from finaletools.frag.wps import wps
from finaletools.frag.agg_wps import aggregate_wps
from finaletools.frag.delfi import delfi
