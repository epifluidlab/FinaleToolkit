"""
Tests for finaletoolkit.utils.frag_generator() and
finaletoolkit.utils.frag_array(), which reads pair-end reads from
several genomic file formats. This function is the backbone of almost
all finaletoolkit features.
"""

import numpy as np
import pytest

from finaletoolkit.utils.utils import frag_generator, frag_array, overlaps

class TestFragGenerator:
    def test_bam(self, request):
        """
        See if frag_generator runs when opening a BAM file and reads the
        right number of reads
        """
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        frag_gen = frag_generator(
            bam, "12", quality_threshold=0, min_length=0, max_length=9999
        )
        frags = [frag for frag in frag_gen]

        chroms = np.array([chrom for chrom, *_ in frags])
        starts = np.array([start for _, start, *_ in frags])
        stops = np.array([stop for _, _, stop, *_ in frags])

        in_region = np.any(overlaps(np.array(['12']), np.array([34442500]),
                               np.array([34446500]), chroms, starts, stops))
        assert in_region, "Some fragments are outside of region"
        
        assert len(frags) == 17, "Incorrect number of frags"

    def test_frag_gz(self, request):
        """
        See if frag_generator runs when opening a frag.gz file and reads the
        right number of reads
        """
        path = request.path.parent / 'data' / '12.3444.b37.frag.gz'
        frag_gen = frag_generator(
            path, "12", quality_threshold=0, min_length=0, max_length=9999
        )
        frags = [frag for frag in frag_gen]

        chroms = np.array([chrom for chrom, *_ in frags])
        starts = np.array([start for _, start, *_ in frags])
        stops = np.array([stop for _, _, stop, *_ in frags])

        in_region = np.any(overlaps(np.array(['12']), np.array([34442500]),
                               np.array([34446500]), chroms, starts, stops))
        assert in_region, "Some fragments are outside of region"
        
        assert len(frags) == 17, "Incorrect number of frags"

    def test_bed_gz(self, request):
        """
        See if frag_generator runs when opening a bed.gz file and reads the
        right number of reads
        """
        path = request.path.parent / 'data' / '12.3444.b37.frag.bed.gz'
        with pytest.warns(UserWarning):
            frag_gen = frag_generator(
                path, "12", quality_threshold=0, min_length=0, max_length=9999
            )
            frags = [frag for frag in frag_gen]

        chroms = np.array([chrom for chrom, *_ in frags])
        starts = np.array([start for _, start, *_ in frags])
        stops = np.array([stop for _, _, stop, *_ in frags])

        in_region = np.any(overlaps(np.array(['12']), np.array([34442500]),
                               np.array([34446500]), chroms, starts, stops))
        assert in_region, "Some fragments are outside of region"
        
        assert len(frags) == 17, "Incorrect number of frags"

    def test_detailed(self,request):
        """
        See if exact fragments are included.
        """
        interval_file = request.path.parent / 'data' / '12.3444.b37.frag.gz'
        contig = "12"
        start=34443119
        stop=34443538
        g = frag_generator(interval_file,contig=contig,start=start,stop=stop)
        expected = [
            ('12', 34443118, 34443284, 60, True),
            ('12', 34443139, 34443300, 60, True),
            ('12', 34443294, 34443491, 60, True),
            ('12', 34443358, 34443538, 60, False),
        ]
        for i, frag in enumerate(g):
            assert frag == expected[i]
        
class TestFragArray:
    def test_bam(self, request):
        """
        See if frag_array runs when opening a BAM file and reads the
        right number of reads
        """
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        frags = frag_array(
            bam, "12", quality_threshold=0, min_length=0, max_length=9999
        )

        starts = frags['start']
        stops = frags['stop']
        strands = frags['strand']

        in_region = np.any(overlaps(np.array(['12']), np.array([34442500]),
                               np.array([34446500]), np.array(['12']), starts, stops))
        assert in_region, "Some fragments are outside of region"
        
        assert len(frags) == 17, "Incorrect number of frags"

    def test_frag_gz(self, request):
        """
        See if frag_array runs when opening a frag.gz file and reads the
        right number of reads
        """
        path = request.path.parent / 'data' / '12.3444.b37.frag.gz'
        frags = frag_array(
            path, "12", quality_threshold=0, min_length=0, max_length=9999
        )

        starts = frags['start']
        stops = frags['stop']
        strands = frags['strand']

        in_region = np.any(
            overlaps(np.array(['12']), np.array([34442500]),
                     np.array([34446500]), np.array(['12']), starts, stops))
        assert in_region, "Some fragments are outside of region"
        
        assert len(frags) == 17, "Incorrect number of frags"
    