"""
Tests for finaletoolkit.utils.frag_generator() and
finaletoolkit.utils.frag_array(), which reads pair-end reads from
several genomic file formats. This function is the backbone of almost
all finaletoolkit features.
"""

import numpy as np

from finaletoolkit.utils.utils import frag_generator, frag_array, overlaps

class TestFragGenerator:
    def test_bam(self, request):
        """
        See if frag_generator runs when opening a BAM file and reads the
        right number of reads
        """
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        frag_gen = frag_generator(
            bam, "12", quality_threshold=0, fraction_low=0, fraction_high=9999
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
            path, "12", quality_threshold=0, fraction_low=0, fraction_high=9999
        )
        frags = [frag for frag in frag_gen]

        chroms = np.array([chrom for chrom, *_ in frags])
        starts = np.array([start for _, start, *_ in frags])
        stops = np.array([stop for _, _, stop, *_ in frags])

        in_region = np.any(overlaps(np.array(['12']), np.array([34442500]),
                               np.array([34446500]), chroms, starts, stops))
        assert in_region, "Some fragments are outside of region"
        
        assert len(frags) == 17, "Incorrect number of frags"

        


class TestFragArray:
    pass
    