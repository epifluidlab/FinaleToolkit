"""
Tests for finaletoolkit.frag.end_motifs module, which calculates
end-motifs and motif diversity scores genomewide and over intervals.
"""

from finaletoolkit.frag.end_motifs import *

class TestGenomewideEndMotifs:
    def test_end_motifs(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        ref_file = request.path.parent / 'data' / '12_34443100-34446700.2bit'
        motifs = end_motifs(
            bam, ref_file, both_strands=True, quality_threshold=0)
