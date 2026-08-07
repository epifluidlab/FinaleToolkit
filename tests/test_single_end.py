"""
Tests for single-end / long-read support.

Single-end libraries (Oxford Nanopore, PacBio) have no mate, so ``TLEN`` is 0
and the pairing flags are unset.  Without ``single_end=True`` every such read
is discarded twice over: once by the pairing filters and again by the
template-length span reconstruction.
"""

import pysam
import pytest
from click.testing import CliRunner

from finaletoolkit.cli.commands import COMMANDS
from finaletoolkit.cli.main_cli import main_cli
from finaletoolkit.io.alignment import AlignmentWrapper
from finaletoolkit.utils import frag_array, frag_generator, low_quality_read_pairs

# Commands that accept long-read input.  delfi and wps are deliberately absent:
# both are defined around nucleosome-length fragments (delfi hard-filters to
# 100-220bp internally) and would produce meaningless output on long reads.
EXPECTED_SINGLE_END_COMMANDS = {
    "coverage",
    "frag-length-bins",
    "frag-length-intervals",
    "cleavage-profile",
    "end-motifs",
    "interval-end-motifs",
    "breakpoint-motifs",
    "interval-breakpoint-motifs",
}

CONTIG = "chr1"
CONTIG_LENGTH = 249250621


def _write_bam(path, reads):
    """Write a coordinate-sorted, indexed BAM from ``(start, cigar, flag, mapq)``."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (start, cigar, flag, mapq) in enumerate(reads):
            seg = pysam.AlignedSegment()
            seg.query_name = f"read_{i}"
            n_query = sum(length for op, length in cigar if op in (0, 1, 4))
            seg.query_sequence = "A" * n_query
            seg.query_qualities = pysam.qualitystring_to_array("I" * n_query)
            seg.flag = flag
            seg.reference_id = 0
            seg.reference_start = start
            seg.mapping_quality = mapq
            seg.cigar = cigar
            seg.template_length = 0  # single-end => TLEN is 0
            out.write(seg)
    pysam.index(str(path))
    return str(path)


@pytest.fixture
def single_end_bam(tmp_path):
    """Three unpaired forward reads spanning 5000, 40000 and 120000 bp."""
    reads = [
        (1_000_000, [(0, 100), (3, 4_900)], 0, 60),
        (1_500_000, [(0, 100), (3, 39_900)], 0, 60),
        (2_000_000, [(0, 100), (3, 119_900)], 0, 60),
    ]
    return _write_bam(tmp_path / "single_end.bam", reads)


class TestSingleEndFragments:
    def test_dropped_without_flag(self, single_end_bam):
        """Regression guard: the pre-fix behaviour, kept as the default."""
        frags = list(frag_generator(single_end_bam, CONTIG, quality_threshold=30))
        assert frags == []

    def test_extracted_with_flag(self, single_end_bam):
        frags = list(
            frag_generator(
                single_end_bam, CONTIG, quality_threshold=30, single_end=True
            )
        )
        assert [(f[1], f[2]) for f in frags] == [
            (1_000_000, 1_005_000),
            (1_500_000, 1_540_000),
            (2_000_000, 2_120_000),
        ]

    def test_span_is_cigar_aware(self, single_end_bam):
        """The span must come from reference_end, not from the stored sequence."""
        lengths = [
            stop - start
            for _, start, stop, _, _ in frag_generator(
                single_end_bam, CONTIG, quality_threshold=30, single_end=True
            )
        ]
        # Each read stores only 100 query bases but spans far more reference.
        assert lengths == [5_000, 40_000, 120_000]

    def test_soft_clipping_excluded(self, tmp_path):
        """Soft-clipped bases are not part of the aligned span."""
        bam = _write_bam(
            tmp_path / "clipped.bam",
            [(500_000, [(4, 50), (0, 1_000), (4, 50)], 0, 60)],
        )
        frags = list(
            frag_generator(bam, CONTIG, quality_threshold=30, single_end=True)
        )
        assert [(f[1], f[2]) for f in frags] == [(500_000, 501_000)]

    def test_frag_array_threads_flag(self, single_end_bam):
        arr = frag_array(single_end_bam, CONTIG, single_end=True)
        assert len(arr) == 3
        assert arr["stop"][0] - arr["start"][0] == 5_000
        assert frag_array(single_end_bam, CONTIG).shape == (0,)

    def test_length_filters_still_apply(self, single_end_bam):
        frags = list(
            frag_generator(
                single_end_bam,
                CONTIG,
                quality_threshold=30,
                min_length=10_000,
                max_length=100_000,
                single_end=True,
            )
        )
        assert [(f[1], f[2]) for f in frags] == [(1_500_000, 1_540_000)]

    def test_reverse_strand_preserved(self, tmp_path):
        bam = _write_bam(
            tmp_path / "reverse.bam",
            [(100_000, [(0, 2_000)], 16, 60)],  # 0x10 => reverse
        )
        frags = list(
            frag_generator(bam, CONTIG, quality_threshold=30, single_end=True)
        )
        assert len(frags) == 1
        assert frags[0][4] is False  # is_forward


class TestSingleEndQualityFiltersStillApply:
    """single_end relaxes only the pairing checks, nothing else."""

    @pytest.mark.parametrize(
        "flag,mapq,reason",
        [
            (0, 10, "below mapping-quality threshold"),
            (4, 60, "unmapped"),
            (256, 60, "secondary"),
            (1024, 60, "duplicate"),
            (2048, 60, "supplementary"),
            (512, 60, "qcfail"),
        ],
    )
    def test_rejected(self, tmp_path, flag, mapq, reason):
        bam = _write_bam(
            tmp_path / f"flag_{flag}_{mapq}.bam",
            [(300_000, [(0, 1_000)], flag, mapq)],
        )
        frags = list(
            frag_generator(bam, CONTIG, quality_threshold=30, single_end=True)
        )
        assert frags == [], f"read should be rejected: {reason}"

    def test_clean_read_accepted(self, tmp_path):
        bam = _write_bam(
            tmp_path / "clean.bam", [(300_000, [(0, 1_000)], 0, 60)]
        )
        frags = list(
            frag_generator(bam, CONTIG, quality_threshold=30, single_end=True)
        )
        assert len(frags) == 1


class TestLowQualityReadPairsSingleEnd:
    def test_unpaired_read_gated_by_flag(self, single_end_bam):
        with pysam.AlignmentFile(single_end_bam) as handle:
            read = next(iter(handle))
        assert low_quality_read_pairs(read) is True
        assert low_quality_read_pairs(read, single_end=True) is False

    def test_mapq_still_enforced_when_single_end(self, tmp_path):
        bam = _write_bam(tmp_path / "lowmapq.bam", [(1_000, [(0, 500)], 0, 5)])
        with pysam.AlignmentFile(bam) as handle:
            read = next(iter(handle))
        assert low_quality_read_pairs(read, min_mapq=30, single_end=True) is True


class TestPairedEndUnchanged:
    """The default path must behave exactly as before."""

    def test_paired_bam_still_reads(self, request):
        bam = request.path.parent / "data" / "12.3444.b37.bam"
        arr = frag_array(bam, "12", min_length=120, max_length=180)
        assert len(arr) > 0

    def test_read1_only_disabled_under_single_end(self, request):
        """single_end forces read1_only off; there is no read1 in such a library."""
        bam = request.path.parent / "data" / "12.3444.b37.bam"
        wrapper = AlignmentWrapper(str(bam), single_end=True)
        assert wrapper.read1_only is False
        wrapper.close()


class TestCliSurface:
    """--single-end is exposed on exactly the long-read-meaningful commands."""

    def test_expected_commands_have_the_flag(self):
        actual = {
            command.name
            for command in COMMANDS
            if any(param.name == "single_end" for param in command.params)
        }
        assert actual == EXPECTED_SINGLE_END_COMMANDS

    @pytest.mark.parametrize("name", sorted(EXPECTED_SINGLE_END_COMMANDS))
    def test_flag_documented_in_help(self, name):
        result = CliRunner().invoke(main_cli, [name, "--help"])
        assert result.exit_code == 0
        assert "--single-end" in result.output

    @pytest.mark.parametrize("name", ["wps", "delfi"])
    def test_nucleosome_commands_reject_the_flag(self, name):
        """delfi and wps must not silently accept a flag they cannot honour."""
        result = CliRunner().invoke(main_cli, [name, "--single-end"])
        assert result.exit_code != 0
        assert "no such option" in result.output.lower()


class TestCliEndToEnd:
    def test_frag_length_bins_needs_the_flag(self, single_end_bam, tmp_path):
        runner = CliRunner()
        without = tmp_path / "without.tsv"
        result = runner.invoke(
            main_cli,
            ["frag-length-bins", single_end_bam, "-o", str(without)],
        )
        assert result.exit_code == 0
        assert not without.exists(), "no fragments should mean no output file"

        with_flag = tmp_path / "with.tsv"
        result = runner.invoke(
            main_cli,
            [
                "frag-length-bins",
                single_end_bam,
                "--single-end",
                "--max-length",
                "200000",
                "--bin-size",
                "1000",
                "-o",
                str(with_flag),
            ],
        )
        assert result.exit_code == 0, result.output
        rows = with_flag.read_text().splitlines()
        assert rows[0] == "min\tmax\tcount"
        assert sum(int(row.split("\t")[2]) for row in rows[1:]) == 3

    def test_coverage_counts_long_reads(self, single_end_bam, tmp_path):
        intervals = tmp_path / "intervals.bed"
        intervals.write_text(f"{CONTIG}\t0\t3000000\tregion\n")
        out = tmp_path / "cov.bed"
        result = CliRunner().invoke(
            main_cli,
            [
                "coverage",
                single_end_bam,
                str(intervals),
                "--single-end",
                "--max-length",
                "200000",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        # Coverage is scaled by --scale-factor, so the value is a float.
        assert float(out.read_text().split("\t")[4]) == 3.0
