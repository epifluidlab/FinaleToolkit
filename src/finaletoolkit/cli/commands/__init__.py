"""
Subcommand builders for the FinaleToolkit CLI.

Each ``_build_<name>`` function registers one subparser (flags, defaults, help,
and an example invocation) and wires it to its implementation via
``set_defaults(module=..., func=...)``.  :func:`build_subparsers` registers all
of them.

The flag scheme is a clean, consistent redesign (see ``cli/_args.py``).  Every
flag's ``dest`` still equals the implementing function's parameter name, so the
Python API is unchanged.
"""
from __future__ import annotations

import argparse

from .._args import (
    add_bool_flag,
    add_input_file,
    add_intersect_policy,
    add_kmer,
    add_max_length,
    add_min_length,
    add_min_mapq,
    add_output,
    add_reference_option,
    add_strand,
    add_threads,
    add_verbose,
)

__all__ = ["build_subparsers"]

_REFSEQ_HELP = (
    "A .2bit or FASTA (.fa, .fasta, .fna) file for the reference genome "
    "sequence used during alignment."
)


def _ex(text: str) -> str:
    """Format an example-invocation epilog."""
    return f"Example: {text}"


def _parser(subparsers, name, description, example):
    """Create a subparser with consistent formatting.

    ``allow_abbrev=False`` disables argparse's fuzzy prefix matching so that,
    e.g., ``--savgol`` is never treated as an ambiguous abbreviation of
    ``--savgol-window-size``. Users type full flag names (as the tests do).
    """
    return subparsers.add_parser(
        name,
        description=description,
        epilog=_ex(example),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
    )


def _build_coverage(subparsers) -> None:
    p = _parser(
        subparsers,
        "coverage",
        "Calculates fragmentation coverage over intervals defined in a BED file "
        "based on alignment data from a BAM/CRAM/Fragment file.",
        "finaletoolkit coverage sample.bam intervals.bed -o cov.bed",
    )
    add_input_file(p)
    add_reference_option(p)
    p.add_argument(
        "interval_file",
        metavar="REGIONS",
        help="Path to a BED file containing intervals to calculate coverage over.",
    )
    add_output(
        p,
        "BED file of coverage values over the input intervals.",
    )
    p.add_argument(
        "-n",
        "--normalize",
        action="store_true",
        help="Normalize output by total coverage (and apply --scale-factor if "
        "given). May increase execution time for high-throughput data.",
    )
    p.add_argument(
        "--scale-factor",
        dest="scale_factor",
        metavar="X",
        default=1.0,
        type=float,
        help="Scale factor for coverage values. Default is 1.",
    )
    add_min_length(p, 0, "Minimum fragment length to include in coverage.")
    add_max_length(p, None, "Maximum fragment length to include in coverage.")
    add_intersect_policy(p)
    add_min_mapq(p, 30)
    add_threads(p)
    add_verbose(p)
    p.set_defaults(module="finaletoolkit.frag._coverage", func="coverage")


def _build_frag_length_bins(subparsers) -> None:
    p = _parser(
        subparsers,
        "frag-length-bins",
        "Retrieves fragment lengths grouped in bins given a BAM/CRAM/Fragment "
        "file.",
        "finaletoolkit frag-length-bins sample.bam --bin-size 5 -o bins.tsv",
    )
    add_input_file(p)
    add_reference_option(p)
    p.add_argument(
        "-c",
        "--contig",
        type=str,
        help="Contig/chromosome to select fragments from. (Required with "
        "--start or --stop.)",
    )
    p.add_argument(
        "-S",
        "--start",
        type=int,
        help="0-based left-most coordinate of the interval to select fragments "
        "from. (Must also specify --contig.)",
    )
    p.add_argument(
        "-E",
        "--stop",
        type=int,
        help="1-based right-most coordinate of the interval to select fragments "
        "from. (Must also specify --contig.)",
    )
    add_min_length(p, 0, "Minimum fragment length to include.")
    add_max_length(p, None, "Maximum fragment length to include.")
    add_intersect_policy(p)
    p.add_argument(
        "--bin-size",
        dest="bin_size",
        metavar="BP",
        type=int,
        default=1,
        help="Width of the bins fragment lengths are grouped into.",
    )
    add_output(p, "TSV of fragment lengths binned by the specified bin size.")
    p.add_argument(
        "--summary-stats",
        dest="summary_stats",
        action="store_true",
        help="Append summary statistics as comment lines (e.g. #max: 100) to "
        "the output TSV.",
    )
    p.add_argument(
        "--short-threshold",
        dest="short_fraction",
        metavar="BP",
        default=None,
        type=int,
        help="If set, include a short fraction (fragments <= this length) in "
        "the summary statistics.",
    )
    p.add_argument(
        "--histogram",
        dest="histogram_path",
        metavar="PNG",
        default=None,
        help="If set, also render a histogram to this PNG path.",
    )
    add_min_mapq(p, 30)
    add_verbose(p)
    p.set_defaults(
        module="finaletoolkit.frag._frag_length", func="frag_length_bins"
    )


def _build_frag_length_intervals(subparsers) -> None:
    p = _parser(
        subparsers,
        "frag-length-intervals",
        "Retrieves fragment length summary statistics over intervals defined in "
        "a BED file based on alignment data from a BAM/CRAM/Fragment file.",
        "finaletoolkit frag-length-intervals sample.bam intervals.bed -o out.bed",
    )
    add_input_file(p)
    add_reference_option(p)
    p.add_argument(
        "interval_file",
        metavar="REGIONS",
        help="Path to a BED file of intervals to summarize fragment lengths over.",
    )
    add_min_length(p, 0, "Minimum fragment length to include.")
    add_max_length(p, None, "Maximum fragment length to include.")
    add_intersect_policy(p)
    add_output(
        p,
        "BED file of fragment-length summary statistics (mean, median, st. dev, "
        "min, max) over the input intervals.",
    )
    p.add_argument(
        "--short-threshold",
        dest="short_reads",
        metavar="BP",
        default=150,
        type=int,
        help="Length cutoff (bp) for the short-read fraction. Default is 150.",
    )
    add_min_mapq(p, 30)
    add_threads(p)
    add_verbose(p)
    p.set_defaults(
        module="finaletoolkit.frag._frag_length", func="frag_length_intervals"
    )


def _build_cleavage_profile(subparsers) -> None:
    p = _parser(
        subparsers,
        "cleavage-profile",
        "Calculates cleavage proportion over intervals defined in a BED file "
        "based on alignment data from a BAM/CRAM/Fragment file.",
        "finaletoolkit cleavage-profile sample.bam cpg.bed hg38.chrom.sizes "
        "--pad-left 5 --pad-right 5 -o cleavage.bw",
    )
    add_input_file(p)
    add_reference_option(p)
    p.add_argument(
        "interval_file",
        metavar="REGIONS",
        help="Path to a BED file of intervals to compute cleavage proportion over.",
    )
    p.add_argument(
        "chrom_sizes",
        metavar="CHROM_SIZES",
        help="A .chrom.sizes file containing chromosome names and sizes.",
    )
    add_output(p, "bigWig file of cleavage proportion over the input intervals.")
    add_min_length(p, 0, "Minimum fragment length to include.")
    add_max_length(p, None, "Maximum fragment length to include.")
    add_min_mapq(p, 20)
    p.add_argument(
        "--pad-left",
        dest="left",
        metavar="BP",
        default=0,
        type=int,
        help="Base pairs to subtract from each start coordinate. Useful for BED "
        "files containing only CpG coordinates. Default is 0.",
    )
    p.add_argument(
        "--pad-right",
        dest="right",
        metavar="BP",
        default=0,
        type=int,
        help="Base pairs to add to each stop coordinate. Useful for BED files "
        "containing only CpG coordinates. Default is 0.",
    )
    add_threads(p)
    add_verbose(p)
    p.set_defaults(
        module="finaletoolkit.frag._cleavage_profile",
        func="multi_cleavage_profile",
    )


def _build_wps(subparsers) -> None:
    p = _parser(
        subparsers,
        "wps",
        "Calculates Windowed Protection Score (WPS) over intervals defined in a "
        "BED file based on alignment data from a BAM/CRAM/Fragment file.",
        "finaletoolkit wps sample.bam tss.bed --chrom-sizes hg38.chrom.sizes "
        "-o wps.bw",
    )
    add_input_file(p)
    add_reference_option(p)
    p.add_argument(
        "site_bed",
        metavar="REGIONS",
        help="Path to a BED file of sites to calculate WPS over, sorted by "
        "contig then start.",
    )
    p.add_argument(
        "--chrom-sizes",
        dest="chrom_sizes",
        metavar="CHROM_SIZES",
        help="A .chrom.sizes file containing chromosome names and sizes.",
    )
    add_output(p, "bigWig file of WPS results over the input intervals.")
    p.add_argument(
        "-i",
        "--interval-size",
        dest="interval_size",
        metavar="BP",
        default=5000,
        type=int,
        help="Size in bp of the windows (centered on each site) to calculate WPS "
        "over. Default is 5000.",
    )
    p.add_argument(
        "-W",
        "--window-size",
        dest="window_size",
        metavar="BP",
        default=120,
        type=int,
        help="Size of the sliding window used to calculate WPS scores. Default "
        "is 120.",
    )
    add_min_length(
        p, 120, "Minimum fragment length to include. Default is 120 (L-WPS)."
    )
    add_max_length(
        p, 180, "Maximum fragment length to include. Default is 180 (L-WPS)."
    )
    add_min_mapq(p, 30)
    add_threads(p)
    add_verbose(p)
    p.set_defaults(module="finaletoolkit.frag", func="multi_wps")


def _build_adjust_wps(subparsers) -> None:
    p = _parser(
        subparsers,
        "adjust-wps",
        "Adjusts raw Windowed Protection Score (WPS) by applying a median filter "
        "and Savitsky-Golay filter.",
        "finaletoolkit adjust-wps wps.bw intervals.bed hg38.chrom.sizes "
        "-o adjusted.bw",
    )
    p.add_argument(
        "input_file",
        metavar="INPUT",
        help="A bigWig file of raw WPS results over the input intervals.",
    )
    p.add_argument(
        "interval_file",
        metavar="REGIONS",
        help="Path to a BED file of intervals WPS was calculated over.",
    )
    p.add_argument(
        "chrom_sizes",
        metavar="CHROM_SIZES",
        help="A .chrom.sizes file containing chromosome names and sizes.",
    )
    add_output(p, "bigWig file of adjusted WPS results over the input intervals.")
    p.add_argument(
        "-i",
        "--interval-size",
        dest="interval_size",
        metavar="BP",
        default=5000,
        type=int,
        help="Size in bp of each interval in the interval file.",
    )
    p.add_argument(
        "-m",
        "--median-window-size",
        dest="median_window_size",
        metavar="BP",
        default=1000,
        type=int,
        help="Size of the median/mean filter window used to adjust WPS scores.",
    )
    p.add_argument(
        "--savgol-window-size",
        dest="savgol_window_size",
        metavar="BP",
        default=21,
        type=int,
        help="Size of the Savitsky-Golay filter window.",
    )
    p.add_argument(
        "--savgol-poly-deg",
        dest="savgol_poly_deg",
        metavar="DEG",
        default=2,
        type=int,
        help="Degree of the Savitsky-Golay filter polynomial.",
    )
    add_bool_flag(
        p,
        "savgol",
        "savgol",
        True,
        "Apply Savitsky-Golay filtering to the adjusted WPS (use --no-savgol to "
        "disable).",
    )
    p.add_argument(
        "--mean",
        action="store_true",
        help="Use a mean filter instead of a median filter.",
    )
    p.add_argument(
        "--subtract-edges",
        dest="subtract_edges",
        action="store_true",
        help="Subtract the mean of each interval's edges (see --edge-size) "
        "before filtering.",
    )
    p.add_argument(
        "--edge-size",
        dest="edge_size",
        metavar="BP",
        default=500,
        type=int,
        help="Edge width subtracted from each interval end when "
        "--subtract-edges is set. Default is 500.",
    )
    add_threads(p)
    add_verbose(p)
    p.set_defaults(module="finaletoolkit.frag._adjust_wps", func="adjust_wps")


def _build_delfi(subparsers) -> None:
    p = _parser(
        subparsers,
        "delfi",
        "Calculates DELFI features over the genome: (GC-corrected) short "
        "fragments, long fragments, DELFI ratio, and total fragments.",
        "finaletoolkit delfi sample.bam autosomes.chrom.sizes hg19.2bit bins.bed "
        "-g hg19 -o delfi.tsv",
    )
    add_input_file(p)
    p.add_argument(
        "chrom_sizes",
        metavar="CHROM_SIZES",
        help="Tab-delimited file of chrom name and integer length. Use only "
        "autosomes to replicate the original scripts.",
    )
    p.add_argument(
        "reference_file",
        metavar="REFERENCE",
        help="A .2bit or FASTA (.fa, .fasta, .fna) reference genome file.",
    )
    p.add_argument(
        "bins_file",
        metavar="BINS",
        help="A BED file of bins over which to calculate DELFI. Use 100kb "
        "autosomal bins to replicate Cristiano et al.",
    )
    p.add_argument(
        "-b",
        "--blacklist",
        dest="blacklist_file",
        metavar="BED",
        help="BED file of regions to ignore when calculating DELFI.",
    )
    p.add_argument(
        "-g",
        "--gap-file",
        dest="gap_file",
        metavar="GAPS",
        help='Telomere/centromere annotations: a BED4 file (4th column '
        '"centromere"/"telomere"/"short arm") or a genome name '
        '(hg19/b37/hg38/GRCh38).',
    )
    add_output(
        p, 'Output file (.bed/.bed.gz/.tsv/.csv). "-" writes to stdout.'
    )
    p.add_argument(
        "--no-gc-correct",
        dest="no_gc_correct",
        action="store_true",
        default=False,
        help="Skip GC correction.",
    )
    add_bool_flag(
        p,
        "remove-nocov",
        "remove_nocov",
        True,
        "Remove the two hg19 regions with no coverage (use --no-remove-nocov "
        "for non-hg19 reference genomes).",
    )
    add_bool_flag(
        p,
        "merge-bins",
        "merge_bins",
        True,
        "Merge input bins to 5Mb (use --no-merge-bins to keep the input bins).",
    )
    p.add_argument(
        "--merge-size",
        dest="window_size",
        metavar="BP",
        default=5000000,
        type=int,
        help="Target size of merged genomic intervals. Default is 5000000.",
    )
    add_min_mapq(p, 30)
    add_threads(p)
    add_verbose(p)
    p.set_defaults(module="finaletoolkit.frag._delfi", func="delfi")


def _build_delfi_gc_correct(subparsers) -> None:
    p = _parser(
        subparsers,
        "delfi-gc-correct",
        "Performs GC-correction on raw DELFI data. Deprecated: the delfi command "
        "performs GC correction by default.",
        "finaletoolkit delfi-gc-correct raw_delfi.bed -o corrected.tsv",
    )
    p.add_argument(
        "input_file",
        metavar="INPUT",
        help='BED file of raw DELFI data with columns "contig", "start", '
        '"stop", "arm", "short", "long", "gc", "num_frags", "ratio".',
    )
    add_output(p, 'BED of GC-corrected DELFI fractions. "-" writes to stdout.')
    p.add_argument(
        "--header-lines",
        dest="header_lines",
        metavar="N",
        default=1,
        type=int,
        help="Number of header lines in the input BED.",
    )
    add_verbose(p)
    p.set_defaults(
        module="finaletoolkit.frag._delfi_gc_correct", func="cli_delfi_gc_correct"
    )


def _build_end_motifs(subparsers) -> None:
    p = _parser(
        subparsers,
        "end-motifs",
        "Measures frequency of k-mer 5' end motifs.",
        "finaletoolkit end-motifs sample.bam hg38.2bit -k 4 -o motifs.tsv",
    )
    add_input_file(p)
    p.add_argument("refseq_file", metavar="REFERENCE", help=_REFSEQ_HELP)
    add_kmer(p, 4)
    add_min_length(p, 50, "Minimum fragment length to include. Default is 50.")
    add_max_length(p, None, "Maximum fragment length to include.")
    add_strand(p, "end")
    add_output(p, 'TSV of k-mer frequencies. "-" writes to stdout.')
    add_min_mapq(p, 20)
    add_threads(p)
    add_verbose(p)
    p.set_defaults(module="finaletoolkit.frag._end_motifs", func="end_motifs")


def _build_interval_end_motifs(subparsers) -> None:
    p = _parser(
        subparsers,
        "interval-end-motifs",
        "Measures frequency of k-mer 5' end motifs in each region of a BED file "
        "and writes a table.",
        "finaletoolkit interval-end-motifs sample.bam hg38.2bit intervals.bed "
        "-o interval_motifs.tsv",
    )
    add_input_file(p)
    p.add_argument("refseq_file", metavar="REFERENCE", help=_REFSEQ_HELP)
    p.add_argument(
        "intervals",
        metavar="REGIONS",
        help="Path to a BED file of intervals to retrieve end-motif frequencies "
        "over.",
    )
    add_kmer(p, 4)
    add_min_length(p, 50, "Minimum fragment length to include. Default is 50.")
    add_max_length(p, None, "Maximum fragment length to include.")
    add_strand(p, "end")
    add_output(p, "Path to a TSV or CSV file to write end-motif frequencies to.")
    add_min_mapq(p, 20)
    add_threads(p)
    add_verbose(p)
    p.set_defaults(
        module="finaletoolkit.frag._end_motifs", func="interval_end_motifs"
    )


def _build_breakpoint_motifs(subparsers) -> None:
    p = _parser(
        subparsers,
        "breakpoint-motifs",
        "Measures frequency of k-mer breakpoint motifs.",
        "finaletoolkit breakpoint-motifs sample.bam hg38.2bit -k 6 -o bp.tsv",
    )
    add_input_file(p)
    p.add_argument("refseq_file", metavar="REFERENCE", help=_REFSEQ_HELP)
    add_kmer(p, 6)
    add_min_length(p, 50, "Minimum fragment length to include. Default is 50.")
    add_max_length(p, None, "Maximum fragment length to include.")
    add_strand(p, "breakpoint")
    add_output(p, 'TSV of k-mer frequencies. "-" writes to stdout.')
    add_min_mapq(p, 20)
    add_threads(p)
    add_verbose(p)
    p.set_defaults(
        module="finaletoolkit.frag._breakpoint_motifs", func="breakpoint_motifs"
    )


def _build_interval_breakpoint_motifs(subparsers) -> None:
    p = _parser(
        subparsers,
        "interval-breakpoint-motifs",
        "Measures frequency of k-mer 5' breakpoint motifs in each region of a "
        "BED file and writes a table.",
        "finaletoolkit interval-breakpoint-motifs sample.bam hg38.2bit "
        "intervals.bed -o interval_bp.tsv",
    )
    add_input_file(p)
    p.add_argument("refseq_file", metavar="REFERENCE", help=_REFSEQ_HELP)
    p.add_argument(
        "intervals",
        metavar="REGIONS",
        help="Path to a BED file of intervals to retrieve breakpoint-motif "
        "frequencies over.",
    )
    add_kmer(p, 6)
    add_min_length(p, 50, "Minimum fragment length to include. Default is 50.")
    add_max_length(p, None, "Maximum fragment length to include.")
    add_strand(p, "breakpoint")
    add_output(
        p, "Path to a TSV or CSV file to write breakpoint-motif frequencies to."
    )
    add_min_mapq(p, 20)
    add_threads(p)
    add_verbose(p)
    p.set_defaults(
        module="finaletoolkit.frag._breakpoint_motifs",
        func="interval_breakpoint_motifs",
    )


def _build_mds(subparsers) -> None:
    p = _parser(
        subparsers,
        "mds",
        "Reads k-mer frequencies from a file and calculates a motif diversity "
        "score (MDS) using normalized Shannon entropy (Jiang et al., 2020).",
        "finaletoolkit mds motifs.tsv",
    )
    p.add_argument(
        "file_path",
        metavar="INPUT",
        nargs="?",
        default="-",
        help="Tabular file with a k-mer column and a frequency column. Reads "
        "from stdin by default.",
    )
    p.add_argument(
        "-s", "--sep", default="\t", help="Field separator. Default is a tab."
    )
    p.add_argument(
        "--header",
        default=0,
        type=int,
        help="Number of header rows to ignore. Default is 0.",
    )
    p.set_defaults(module="finaletoolkit.frag._end_motifs", func="_cli_mds")


def _build_interval_mds(subparsers) -> None:
    p = _parser(
        subparsers,
        "interval-mds",
        "Reads interval k-mer frequencies and calculates a motif diversity "
        "score (MDS) per interval (Jiang et al., 2020).",
        "finaletoolkit interval-mds interval_motifs.tsv mds.bed",
    )
    p.add_argument(
        "file_path",
        metavar="INPUT",
        nargs="?",
        default="-",
        help="Tabular file of interval k-mer frequencies. Reads from stdin by "
        "default.",
    )
    p.add_argument(
        "-s", "--sep", default="\t", help="Field separator. Default is a tab."
    )
    p.add_argument(
        "file_out",
        metavar="OUTPUT",
        default="-",
        help="Output BED/BEDGraph file of per-interval MDS.",
    )
    p.add_argument(
        "--header",
        default=0,
        type=int,
        help="Number of header rows to ignore. Default is 0.",
    )
    p.set_defaults(
        module="finaletoolkit.frag._end_motifs", func="_cli_interval_mds"
    )


def _build_filter_file(subparsers) -> None:
    p = _parser(
        subparsers,
        "filter-file",
        "Filters a BED/BAM/CRAM file so reads/intervals are properly paired, "
        "exceed a MAPQ, are read1, are not secondary/supplementary/qcfail, are "
        "within/excluding given regions, and share a reference with their mate.",
        "finaletoolkit filter-file sample.bam -q 30 --min-length 100 "
        "-o filtered.bam",
    )
    p.add_argument("input_file", metavar="INPUT", help="Path to a BAM/CRAM/BED file.")
    p.add_argument(
        "-w",
        "--whitelist",
        dest="whitelist_file",
        metavar="BED",
        default=None,
        help="Only keep alignments overlapping intervals in this BED file.",
    )
    p.add_argument(
        "-b",
        "--blacklist",
        dest="blacklist_file",
        metavar="BED",
        default=None,
        help="Only keep alignments outside intervals in this BED file.",
    )
    add_output(p, "Output BED/BAM/CRAM file path.")
    add_min_mapq(p, 30)
    add_min_length(p, None, "Minimum fragment length to include.")
    add_max_length(p, None, "Maximum fragment length to include.")
    add_intersect_policy(p)
    add_threads(p)
    add_verbose(p)
    p.set_defaults(module="finaletoolkit.utils", func="filter_file")


def _build_agg_bw(subparsers) -> None:
    p = _parser(
        subparsers,
        "agg-bw",
        "Aggregates a bigWig signal over constant-length intervals defined in a "
        "BED file.",
        "finaletoolkit agg-bw wps.bw intervals.bed -m 120 -o agg.wig",
    )
    p.add_argument(
        "input_file",
        metavar="INPUT",
        help="A bigWig file of signal over the input intervals.",
    )
    p.add_argument(
        "interval_file",
        metavar="REGIONS",
        help="Path to a BED file of intervals the signal was calculated over.",
    )
    add_output(p, "Wiggle file of the aggregate signal over the input intervals.")
    p.add_argument(
        "-m",
        "--median-window-size",
        dest="median_window_size",
        metavar="BP",
        default=1,
        type=int,
        help="Median filter window used to aggregate scores. Set to 120 when "
        "aggregating WPS signals.",
    )
    p.add_argument(
        "--mean",
        action="store_true",
        help="Use the mean instead of the sum.",
    )
    add_verbose(p)
    p.set_defaults(module="finaletoolkit.utils._agg_bw", func="agg_bw")


def _build_gap_bed(subparsers) -> None:
    p = subparsers.add_parser(
        "gap-bed",
        description="Creates a BED4 file of centromeres, telomeres, and "
        "short-arm intervals (the UCSC hg19 'gaps' track; Kent et al., 2002). "
        "Supports hg19, b37, human_g1k_v37, hg38, and GRCh38.",
        epilog="Gap is used liberally; for hg38/GRCh38 it may refer to regions "
        "that no longer have gaps in the reference.\n"
        + _ex("finaletoolkit gap-bed hg19 hg19_gaps.bed"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
    )
    p.add_argument(
        "reference_genome",
        metavar="GENOME",
        choices=["hg19", "b37", "human_g1k_v37", "hg38", "GRCh38"],
        help="Reference genome to provide gaps for.",
    )
    p.add_argument(
        "output_file",
        metavar="OUTPUT",
        help='Path to write the BED file to. "-" writes to stdout.',
    )
    p.set_defaults(module="finaletoolkit.genome.gaps", func="_cli_gap_bed")


# Ordered builders (matches the original subcommand registration order).
_BUILDERS = [
    _build_coverage,
    _build_frag_length_bins,
    _build_frag_length_intervals,
    _build_cleavage_profile,
    _build_wps,
    _build_adjust_wps,
    _build_delfi,
    _build_delfi_gc_correct,
    _build_end_motifs,
    _build_interval_end_motifs,
    _build_breakpoint_motifs,
    _build_interval_breakpoint_motifs,
    _build_mds,
    _build_interval_mds,
    _build_filter_file,
    _build_agg_bw,
    _build_gap_bed,
]


def build_subparsers(subparsers) -> None:
    """Register every FinaleToolkit subcommand on ``subparsers``."""
    for builder in _BUILDERS:
        builder(subparsers)
