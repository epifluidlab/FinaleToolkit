"""
Tests for finaletoolkit.genome.gaps
"""

import os

import pytest

from finaletoolkit.genome.gaps import *

class TestGapBedCLI:
    def test_hg19(self, tmp_path):
        dest = tmp_path / 'hg19.gaps.bed'

        exit_status = os.system(f'finaletoolkit gap-bed hg19 {dest}')
        assert exit_status == 0

        file_exists = os.path.isfile(dest)
        assert file_exists
        

    def test_b37(self, tmp_path):
        dest = tmp_path / 'b37.gaps.bed'

        exit_status = os.system(f'finaletoolkit gap-bed b37 {dest}')
        assert exit_status == 0

        file_exists = os.path.isfile(dest)
        assert file_exists
        

    def test_human_g1k_v37(self, tmp_path):
        dest = tmp_path / 'human_g1k_v37.gaps.bed'

        exit_status = os.system(f'finaletoolkit gap-bed human_g1k_v37 {dest}')
        assert exit_status == 0

        file_exists = os.path.isfile(dest)
        assert file_exists
        

    def test_hg38(self, tmp_path):
        dest = tmp_path / 'hg38.gaps.bed'

        exit_status = os.system(f'finaletoolkit gap-bed hg38 {dest}')
        assert exit_status == 0

        file_exists = os.path.isfile(dest)
        assert file_exists
        

    def test_GRCh38(self, tmp_path):
        dest = tmp_path / 'GRCh38.gaps.bed'

        exit_status = os.system(f'finaletoolkit gap-bed GRCh38 {dest}')
        assert exit_status == 0

        file_exists = os.path.isfile(dest)
        assert file_exists

class TestGapBed:
    def test_hg19(self, tmp_path):
        dest = tmp_path / 'hg19.gaps.bed'

        ucsc_hg19_gap_bed(dest)

        file_exists = os.path.isfile(dest)
        assert file_exists

    def test_b37(self, tmp_path):
        dest = tmp_path / 'b37.gaps.bed'

        b37_gap_bed(dest)

        file_exists = os.path.isfile(dest)
        assert file_exists

    def test_hg38(self, tmp_path):
        dest = tmp_path / 'hg19.gaps.bed'

        ucsc_hg38_gap_bed(dest)

        file_exists = os.path.isfile(dest)
        assert file_exists

class TestGenomeGaps:
    pass
