"""
Tests for finaletoolkit.utils
"""

import os
import filecmp
import difflib
import pysam
from numpy import array
from numpy.testing import assert_array_equal    

import pytest


from finaletoolkit.utils import *

class TestUtils:
    def test_chrom_sizes_to_dict(self, request):
        chrom_sizes = request.path.parent / 'data' / 'b37.chrom.sizes'
        assert chrom_sizes_to_dict(chrom_sizes) == {'1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276, '5': 180915260, '6': 171115067, '7': 159138663, '8': 146364022, '9': 141213431, '10': 135534747, '11': 135006516, '12': 133851895, '13': 115169878, '14': 107349540, '15': 102531392, '16': 90354753, '17': 81195210, '18': 78077248, '19': 59128983, '20': 63025520, '21': 48129895, '22': 51304566, 'X': 155270560, 'Y': 59373566, 'MT': 16569, 'GL000207.1': 4262, 'GL000226.1': 15008, 'GL000229.1': 19913, 'GL000231.1': 27386, 'GL000210.1': 27682, 'GL000239.1': 33824, 'GL000235.1': 34474, 'GL000201.1': 36148, 'GL000247.1': 36422, 'GL000245.1': 36651, 'GL000197.1': 37175, 'GL000203.1': 37498, 'GL000246.1': 38154, 'GL000249.1': 38502, 'GL000196.1': 38914, 'GL000248.1': 39786, 'GL000244.1': 39929, 'GL000238.1': 39939, 'GL000202.1': 40103, 'GL000234.1': 40531, 'GL000232.1': 40652, 'GL000206.1': 41001, 'GL000240.1': 41933, 'GL000236.1': 41934, 'GL000241.1': 42152, 'GL000243.1': 43341, 'GL000242.1': 43523, 'GL000230.1': 43691, 'GL000237.1': 45867, 'GL000233.1': 45941, 'GL000204.1': 81310, 'GL000198.1': 90085, 'GL000208.1': 92689, 'GL000191.1': 106433, 'GL000227.1': 128374, 'GL000228.1': 129120, 'GL000214.1': 137718, 'GL000221.1': 155397, 'GL000209.1': 159169, 'GL000218.1': 161147, 'GL000220.1': 161802, 'GL000213.1': 164239, 'GL000211.1': 166566, 'GL000199.1': 169874, 'GL000217.1': 172149, 'GL000216.1': 172294, 'GL000215.1': 172545, 'GL000205.1': 174588, 'GL000219.1': 179198, 'GL000224.1': 179693, 'GL000223.1': 180455, 'GL000195.1': 182896, 'GL000212.1': 186858, 'GL000222.1': 186861, 'GL000200.1': 187035, 'GL000193.1': 189789, 'GL000194.1': 191469, 'GL000225.1': 211173, 'GL000192.1': 547496, 'NC_007605': 171823}

    def test_chrom_sizes_to_list(self, request):
        chrom_sizes = request.path.parent / 'data' / 'b37.chrom.sizes'
        assert chrom_sizes_to_list(chrom_sizes) == [('1', 249250621), ('2', 243199373), ('3', 198022430), ('4', 191154276), ('5', 180915260), ('6', 171115067), ('7', 159138663), ('8', 146364022), ('9', 141213431), ('10', 135534747), ('11', 135006516), ('12', 133851895), ('13', 115169878), ('14', 107349540), ('15', 102531392), ('16', 90354753), ('17', 81195210), ('18', 78077248), ('19', 59128983), ('20', 63025520), ('21', 48129895), ('22', 51304566), ('X', 155270560), ('Y', 59373566), ('MT', 16569), ('GL000207.1', 4262), ('GL000226.1', 15008), ('GL000229.1', 19913), ('GL000231.1', 27386), ('GL000210.1', 27682), ('GL000239.1', 33824), ('GL000235.1', 34474), ('GL000201.1', 36148), ('GL000247.1', 36422), ('GL000245.1', 36651), ('GL000197.1', 37175), ('GL000203.1', 37498), ('GL000246.1', 38154), ('GL000249.1', 38502), ('GL000196.1', 38914), ('GL000248.1', 39786), ('GL000244.1', 39929), ('GL000238.1', 39939), ('GL000202.1', 40103), ('GL000234.1', 40531), ('GL000232.1', 40652), ('GL000206.1', 41001), ('GL000240.1', 41933), ('GL000236.1', 41934), ('GL000241.1', 42152), ('GL000243.1', 43341), ('GL000242.1', 43523), ('GL000230.1', 43691), ('GL000237.1', 45867), ('GL000233.1', 45941), ('GL000204.1', 81310), ('GL000198.1', 90085), ('GL000208.1', 92689), ('GL000191.1', 106433), ('GL000227.1', 128374), ('GL000228.1', 129120), ('GL000214.1', 137718), ('GL000221.1', 155397), ('GL000209.1', 159169), ('GL000218.1', 161147), ('GL000220.1', 161802), ('GL000213.1', 164239), ('GL000211.1', 166566), ('GL000199.1', 169874), ('GL000217.1', 172149), ('GL000216.1', 172294), ('GL000215.1', 172545), ('GL000205.1', 174588), ('GL000219.1', 179198), ('GL000224.1', 179693), ('GL000223.1', 180455), ('GL000195.1', 182896), ('GL000212.1', 186858), ('GL000222.1', 186861), ('GL000200.1', 187035), ('GL000193.1', 189789), ('GL000194.1', 191469), ('GL000225.1', 211173), ('GL000192.1', 547496), ('NC_007605', 171823)]

    def test_frag_array(self, request):
        interval_file = request.path.parent / 'data' / '12.3444.b37.frag.gz'
        contig = "12"
        arr = frag_array(interval_file, contig, min_length=120, max_length=180)
        assert_array_equal(
            arr,
            array([
                (34443118, 34443284,  True), (34443139, 34443300,  True),
                (34443358, 34443538, False), (34443483, 34443660,  True),
                (34444089, 34444252,  True), (34444696, 34444863,  True),
                (34444954, 34445075,  True), (34444968, 34445105,  True),
                (34445136, 34445288,  True), (34445511, 34445672, False),
                (34445705, 34445852,  True), (34445723, 34445893,  True),
                (34446126, 34446261, False), (34446486, 34446653,  True)],
                dtype=[('start', '<i8'), ('stop', '<i8'), ('strand', '?')]
              )
        )        
    
    def test_frags_in_region(self, request):
        interval_file = str(
              request.path.parent / 'data' / '12.3444.b37.frag.gz') # Numba can't work with PosixPath        
        start = 34443119
        stop = 34445075
        contig = "12"
        arr = frag_array(interval_file, contig, min_length=120, max_length=180)
        frags = frags_in_region(arr, start, stop)
        assert_array_equal(
            frags,
            array([
                (34443118, 34443284,  True), (34443139, 34443300,  True),
                (34443358, 34443538, False), (34443483, 34443660,  True),
                (34444089, 34444252,  True), (34444696, 34444863,  True),
                (34444954, 34445075,  True), (34444968, 34445105,  True)],
                dtype=[('start', '<i8'), ('stop', '<i8'), ('strand', '?')]
            )
        )

    def test_low_quality_read_pairs(self, request):
        interval_file = request.path.parent / 'data' / '12.3444.b37.bam'
        read = pysam.AlignmentFile(interval_file)
        assert low_quality_read_pairs(next(read))
        assert low_quality_read_pairs(next(read))
        assert low_quality_read_pairs(next(read))
        assert not low_quality_read_pairs(next(read))


class TestNoneToleranceComparisons:
    """_none_leq/_none_geq/_none_eq: a None operand means "unbounded"."""

    def test_none_leq(self):
        assert _none_leq(1, 2)
        assert not _none_leq(2, 1)
        assert _none_leq(None, 2)
        assert _none_leq(1, None)
        assert _none_leq(None, None)

    def test_none_geq(self):
        assert _none_geq(2, 1)
        assert not _none_geq(1, 2)
        assert _none_geq(None, 2)
        assert _none_geq(1, None)
        assert _none_geq(None, None)

    def test_none_eq(self):
        assert _none_eq(3, 3)
        assert not _none_eq(3, 4)
        assert _none_eq(None, 3)
        assert _none_eq(3, None)
        assert _none_eq(None, None)


class TestMergeIntervals:
    """_merge_overlapping_intervals and friends, used to collapse a BED file's
    intervals per contig before further processing."""

    def test_merge_overlapping_intervals(self):
        intervals = [(10, 20), (15, 25), (30, 40), (100, 200)]
        assert _merge_overlapping_intervals(intervals) == [
            (10, 25), (30, 40), (100, 200)
        ]

    def test_merge_overlapping_intervals_no_overlap(self):
        intervals = [(30, 40), (10, 20)]
        assert _merge_overlapping_intervals(intervals) == [(10, 20), (30, 40)]

    def test_merge_overlapping_intervals_containment(self):
        # A fully-contained interval shouldn't shrink the merged span.
        assert _merge_overlapping_intervals([(10, 100), (20, 30)]) == [(10, 100)]

    def test_merge_overlapping_intervals_empty(self):
        assert _merge_overlapping_intervals([]) == []

    def test_reduce_overlaps_in_file(self, tmp_path):
        bed = tmp_path / "intervals.bed"
        bed.write_text(
            "1\t10\t20\n"
            "1\t15\t25\n"
            "2\t5\t8\n"
        )
        assert _reduce_overlaps_in_file(str(bed)) == {
            "1": [(10, 25)],
            "2": [(5, 8)],
        }

    def test_convert_to_list(self):
        reduced = {"1": [(10, 20), (30, 40)]}
        assert _convert_to_list(reduced) == {
            "1": [["1", 10, 20], ["1", 30, 40]]
        }

    def test_merge_all_intervals(self):
        converted = {
            "1": [["1", 10, 20]],
            "2": [["2", 5, 8], ["2", 50, 60]],
        }
        result = _merge_all_intervals(converted)
        assert result == [["1", 10, 20], ["2", 5, 8], ["2", 50, 60]]

    def test_merge_all_intervals_empty(self):
        assert _merge_all_intervals({}) == []