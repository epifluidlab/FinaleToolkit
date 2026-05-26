"""
Tests for finaletoolkit.frag.wps, finaletoolkit.frag.adjust_wps, and
finaletoolkit.frag.multi_wps
"""

import os
import filecmp
import difflib

import pytest
import numpy as np
import pysam
import pyBigWig as pbw

from finaletoolkit.frag import wps, multi_wps

class TestWPS:
    def test_lwps(self, request):
        bam = request.path.parent / 'data' / '12.3444.b37.bam'
        results = wps(bam, '12', 34444145, 34444155, 133851895, quality_threshold=0)

        assert np.all(results['contig'] == '12')
        assert np.all(
            results['start'] == np.arange(34444145, 34444155)), str(results)
        assert np.all(
            results['wps'] == [-1,-1,-1,-1,-1,1,1,1,1,1]), str(results)


def _make_two_chrom_bam(path: str) -> str:
    """
    Create a coordinate-sorted, indexed BAM with one proper fragment on each
    of two chromosomes.  Header order is ["2", "10"] (numeric), so an
    alphabetically-sorted BED file lists "10" before "2" — the ordering that
    previously caused silent chromosome dropout in multi_wps.

    Fragment layout (read_len=150, frag_len=160, window_size=120):
      chr  start      stop       WPS non-zero range
      2    1_000_000  1_000_160  [1_000_060, 1_000_100)
      10   1_000_000  1_000_160  [1_000_060, 1_000_100)
    """
    chroms = [("2", 100_000_000), ("10", 100_000_000)]
    hdr = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": l} for c, l in chroms],
    })
    read_len = 150
    frag_len = 160  # within default min_length=120, max_length=180
    chrom_idx = {c: i for i, (c, _) in enumerate(chroms)}

    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=hdr) as bam:
        for fid, (chrom, pos) in enumerate([("2", 1_000_000), ("10", 1_000_000)]):
            rid = chrom_idx[chrom]
            r2_pos = pos + frag_len - read_len  # = pos + 10

            r1 = pysam.AlignedSegment(hdr)
            r1.query_name = f"f{fid}"
            r1.reference_id = rid
            r1.reference_start = pos
            r1.cigarstring = f"{read_len}M"
            r1.flag = 99   # paired, proper pair, mate reverse strand, read1
            r1.mapping_quality = 60
            r1.query_sequence = "A" * read_len
            r1.query_qualities = pysam.qualitystring_to_array("I" * read_len)
            r1.next_reference_id = rid
            r1.next_reference_start = r2_pos
            r1.template_length = frag_len
            bam.write(r1)

            r2 = pysam.AlignedSegment(hdr)
            r2.query_name = f"f{fid}"
            r2.reference_id = rid
            r2.reference_start = r2_pos
            r2.cigarstring = f"{read_len}M"
            r2.flag = 147  # paired, proper pair, read reverse strand, read2
            r2.mapping_quality = 60
            r2.query_sequence = "A" * read_len
            r2.query_qualities = pysam.qualitystring_to_array("I" * read_len)
            r2.next_reference_id = rid
            r2.next_reference_start = pos
            r2.template_length = -frag_len
            bam.write(r2)

    pysam.sort("-o", path, unsorted)
    pysam.index(path)
    os.unlink(unsorted)
    return path


class TestMultiWpsBedOrder:
    """
    Regression test for silent chromosome dropout in multi_wps.

    When the input BED is sorted alphabetically ("10" before "2") but the BAM
    header is in numeric order ("2" before "10"), pyBigWig raises RuntimeError
    on out-of-order writes.  The RuntimeError handler previously swallowed the
    error with `continue`, silently dropping all data for the affected
    chromosome.  The fix sorts intervals by header chromosome order before
    passing them to the pool.
    """

    @pytest.fixture(scope="class")
    def two_chrom_bam(self, tmp_path_factory):
        tmp = tmp_path_factory.mktemp("bam_order")
        return _make_two_chrom_bam(str(tmp / "test.bam"))

    def test_all_chroms_present_when_bed_alphabetical(self, two_chrom_bam, tmp_path):
        # BED in alphabetical order: "10" sorts before "2" as a string.
        # Without the interval-sort fix, chromosome "2" entries would be
        # silently dropped because "2" precedes "10" in the BAM header.
        bed = tmp_path / "sites.bed"
        bed.write_text("10\t999500\t1000500\n2\t999500\t1000500\n")
        output = str(tmp_path / "out.bw")

        multi_wps(
            two_chrom_bam,
            str(bed),
            output_file=output,
            interval_size=1000,
            min_length=120,
            max_length=180,
            quality_threshold=0,
        )

        with pbw.open(output, "r") as bw:
            chr2_max = bw.stats("2", 999_000, 1_001_000, type="max")
            chr10_max = bw.stats("10", 999_000, 1_001_000, type="max")

        assert chr2_max[0] is not None, (
            "Chromosome '2' has no BigWig entries — likely dropped due to "
            "out-of-order BED/header mismatch"
        )
        assert chr10_max[0] is not None, (
            "Chromosome '10' has no BigWig entries"
        )