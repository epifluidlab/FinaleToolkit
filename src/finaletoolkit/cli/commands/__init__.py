"""
Click subcommands for the FinaleToolkit CLI.

Each command collects its options into ``**params`` and hands off to
:func:`finaletoolkit.cli._dispatch.run`, which lazily imports the implementing
function (looked up in :data:`COMMAND_TARGETS`) and calls it.  Keeping a single
``name -> (module, func)`` registry means the dispatch target lives in exactly
one place and stays trivially introspectable (the test suite walks it).

The flag scheme is the consistent redesign documented in
:mod:`finaletoolkit.cli._args`; every option's parameter name still equals the
implementing function's argument, so the Python API is unchanged.
"""
from __future__ import annotations

import rich_click as click

from .._args import (
    bool_flag,
    input_argument,
    intersect_policy_option,
    kmer_option,
    max_length_option,
    min_length_option,
    min_mapq_option,
    output_option,
    reference_option,
    strand_option,
    threads_option,
    verbose_option,
)
from .._dispatch import run

__all__ = ["COMMANDS", "COMMAND_TARGETS", "register_commands"]

# Single source of truth: subcommand name -> (implementing module, function).
COMMAND_TARGETS: dict[str, tuple[str, str]] = {
    "coverage": ("finaletoolkit.frag._coverage", "coverage"),
    "frag-length-bins": ("finaletoolkit.frag._frag_length", "frag_length_bins"),
    "frag-length-intervals": (
        "finaletoolkit.frag._frag_length",
        "frag_length_intervals",
    ),
    "cleavage-profile": (
        "finaletoolkit.frag._cleavage_profile",
        "multi_cleavage_profile",
    ),
    "wps": ("finaletoolkit.frag", "multi_wps"),
    "adjust-wps": ("finaletoolkit.frag._adjust_wps", "adjust_wps"),
    "delfi": ("finaletoolkit.frag._delfi", "delfi"),
    "end-motifs": ("finaletoolkit.frag._end_motifs", "end_motifs"),
    "interval-end-motifs": (
        "finaletoolkit.frag._end_motifs",
        "interval_end_motifs",
    ),
    "breakpoint-motifs": (
        "finaletoolkit.frag._breakpoint_motifs",
        "breakpoint_motifs",
    ),
    "interval-breakpoint-motifs": (
        "finaletoolkit.frag._breakpoint_motifs",
        "interval_breakpoint_motifs",
    ),
    "mds": ("finaletoolkit.frag._end_motifs", "_cli_mds"),
    "regional-mds": ("finaletoolkit.frag._end_motifs", "_cli_regional_mds"),
    "filter-file": ("finaletoolkit.utils", "filter_file"),
    "agg-bw": ("finaletoolkit.utils._agg_bw", "agg_bw"),
    "gap-bed": ("finaletoolkit.genome.gaps", "_cli_gap_bed"),
}

_REFSEQ_HELP = (
    "A .2bit or FASTA (.fa, .fasta, .fna) reference genome file used during "
    "alignment."
)


def _dispatch(**params) -> None:
    """Look up the current command's target and run it with ``params``."""
    name = click.get_current_context().command.name
    module, func = COMMAND_TARGETS[name]
    run(module, func, params)


def _ex(text: str) -> str:
    """Format an example-invocation epilog."""
    return f"Example: {text}"


# --------------------------------------------------------------------------- #
# Coverage & fragment length
# --------------------------------------------------------------------------- #
@click.command("coverage", epilog=_ex("finaletoolkit coverage sample.bam intervals.bed -o cov.bed"))
@input_argument()
@click.argument("interval_file", metavar="REGIONS")
@reference_option()
@output_option("BED file of coverage values over the input intervals.")
@click.option(
    "-n",
    "--normalize",
    "normalize",
    is_flag=True,
    help="Normalize output by total coverage (and apply --scale-factor if "
    "given). May increase execution time for high-throughput data.",
)
@click.option(
    "--scale-factor",
    "scale_factor",
    metavar="X",
    default=1.0,
    show_default=True,
    type=float,
    help="Scale factor for coverage values.",
)
@min_length_option(0, "Minimum fragment length to include in coverage.")
@max_length_option(None, "Maximum fragment length to include in coverage.")
@intersect_policy_option()
@min_mapq_option(30)
@threads_option()
@verbose_option()
def coverage(**params):
    """Fragmentation coverage over BED intervals from a BAM/CRAM/Fragment file.

    \b
    INPUT    BAM/CRAM/Fragment file of fragment data.
    REGIONS  BED file of intervals to calculate coverage over.
    """
    _dispatch(**params)


@click.command("frag-length-bins", epilog=_ex("finaletoolkit frag-length-bins sample.bam --bin-size 5 -o bins.tsv"))
@input_argument()
@reference_option()
@click.option(
    "-c",
    "--contig",
    "contig",
    type=str,
    help="Contig/chromosome to select fragments from. (Required with --start "
    "or --stop.)",
)
@click.option(
    "-S",
    "--start",
    "start",
    type=int,
    help="0-based left-most coordinate of the interval. (Requires --contig.)",
)
@click.option(
    "-E",
    "--stop",
    "stop",
    type=int,
    help="1-based right-most coordinate of the interval. (Requires --contig.)",
)
@min_length_option(0, "Minimum fragment length to include.")
@max_length_option(None, "Maximum fragment length to include.")
@intersect_policy_option()
@click.option(
    "--bin-size",
    "bin_size",
    metavar="BP",
    type=int,
    default=1,
    show_default=True,
    help="Width of the bins fragment lengths are grouped into.",
)
@output_option("TSV of fragment lengths binned by the specified bin size.")
@click.option(
    "--summary-stats",
    "summary_stats",
    is_flag=True,
    help="Append summary statistics as comment lines (e.g. #max: 100) to the "
    "output TSV.",
)
@click.option(
    "--short-threshold",
    "short_fraction",
    metavar="BP",
    default=None,
    type=int,
    help="If set, include a short fraction (fragments <= this length) in the "
    "summary statistics.",
)
@click.option(
    "--histogram",
    "histogram_path",
    metavar="PNG",
    default=None,
    help="If set, also render a histogram to this PNG path.",
)
@min_mapq_option(30)
@verbose_option()
def frag_length_bins(**params):
    """Fragment lengths grouped into bins from a BAM/CRAM/Fragment file.

    \b
    INPUT  BAM/CRAM/Fragment file of fragment data.
    """
    _dispatch(**params)


@click.command("frag-length-intervals", epilog=_ex("finaletoolkit frag-length-intervals sample.bam intervals.bed -o out.bed"))
@input_argument()
@click.argument("interval_file", metavar="REGIONS")
@reference_option()
@min_length_option(0, "Minimum fragment length to include.")
@max_length_option(None, "Maximum fragment length to include.")
@intersect_policy_option()
@output_option(
    "BED of fragment-length summary statistics (mean, median, st. dev, min, "
    "max) over the input intervals."
)
@click.option(
    "--short-threshold",
    "short_reads",
    metavar="BP",
    default=150,
    show_default=True,
    type=int,
    help="Length cutoff (bp) for the short-read fraction.",
)
@min_mapq_option(30)
@threads_option()
@verbose_option()
def frag_length_intervals(**params):
    """Fragment-length summary statistics over BED intervals.

    \b
    INPUT    BAM/CRAM/Fragment file of fragment data.
    REGIONS  BED file of intervals to summarize fragment lengths over.
    """
    _dispatch(**params)


# --------------------------------------------------------------------------- #
# WPS & cleavage
# --------------------------------------------------------------------------- #
@click.command("cleavage-profile", epilog=_ex("finaletoolkit cleavage-profile sample.bam cpg.bed hg38.chrom.sizes --pad-left 5 --pad-right 5 -o cleavage.bw"))
@input_argument()
@click.argument("interval_file", metavar="REGIONS")
@click.argument("chrom_sizes", metavar="CHROM_SIZES")
@reference_option()
@output_option("bigWig file of cleavage proportion over the input intervals.")
@min_length_option(0, "Minimum fragment length to include.")
@max_length_option(None, "Maximum fragment length to include.")
@min_mapq_option(20)
@click.option(
    "--pad-left",
    "left",
    metavar="BP",
    default=0,
    show_default=True,
    type=int,
    help="Base pairs to subtract from each start coordinate. Useful for BED "
    "files containing only CpG coordinates.",
)
@click.option(
    "--pad-right",
    "right",
    metavar="BP",
    default=0,
    show_default=True,
    type=int,
    help="Base pairs to add to each stop coordinate. Useful for BED files "
    "containing only CpG coordinates.",
)
@threads_option()
@verbose_option()
def cleavage_profile(**params):
    """Cleavage proportion over BED intervals from a BAM/CRAM/Fragment file.

    \b
    INPUT        BAM/CRAM/Fragment file of fragment data.
    REGIONS      BED file of intervals to compute cleavage proportion over.
    CHROM_SIZES  .chrom.sizes file of chromosome names and sizes.
    """
    _dispatch(**params)


@click.command("wps", epilog=_ex("finaletoolkit wps sample.bam tss.bed --chrom-sizes hg38.chrom.sizes -o wps.bw"))
@input_argument()
@click.argument("site_bed", metavar="REGIONS")
@reference_option()
@click.option(
    "--chrom-sizes",
    "chrom_sizes",
    metavar="CHROM_SIZES",
    help="A .chrom.sizes file containing chromosome names and sizes.",
)
@output_option("bigWig file of WPS results over the input intervals.")
@click.option(
    "-i",
    "--interval-size",
    "interval_size",
    metavar="BP",
    default=5000,
    show_default=True,
    type=int,
    help="Size in bp of the windows (centered on each site) to calculate WPS "
    "over.",
)
@click.option(
    "-W",
    "--window-size",
    "window_size",
    metavar="BP",
    default=120,
    show_default=True,
    type=int,
    help="Size of the sliding window used to calculate WPS scores.",
)
@min_length_option(120, "Minimum fragment length to include (L-WPS).")
@max_length_option(180, "Maximum fragment length to include (L-WPS).")
@min_mapq_option(30)
@threads_option()
@verbose_option()
def wps(**params):
    """Windowed Protection Score (WPS) over BED sites.

    \b
    INPUT    BAM/CRAM/Fragment file of fragment data.
    REGIONS  BED file of sites, sorted by contig then start.
    """
    _dispatch(**params)


@click.command("adjust-wps", epilog=_ex("finaletoolkit adjust-wps wps.bw intervals.bed hg38.chrom.sizes -o adjusted.bw"))
@click.argument("input_file", metavar="INPUT")
@click.argument("interval_file", metavar="REGIONS")
@click.argument("chrom_sizes", metavar="CHROM_SIZES")
@output_option("bigWig file of adjusted WPS results over the input intervals.")
@click.option(
    "-i",
    "--interval-size",
    "interval_size",
    metavar="BP",
    default=5000,
    show_default=True,
    type=int,
    help="Size in bp of each interval in the interval file.",
)
@click.option(
    "-m",
    "--median-window-size",
    "median_window_size",
    metavar="BP",
    default=1000,
    show_default=True,
    type=int,
    help="Size of the median/mean filter window used to adjust WPS scores.",
)
@click.option(
    "--savgol-window-size",
    "savgol_window_size",
    metavar="BP",
    default=21,
    show_default=True,
    type=int,
    help="Size of the Savitsky-Golay filter window.",
)
@click.option(
    "--savgol-poly-deg",
    "savgol_poly_deg",
    metavar="DEG",
    default=2,
    show_default=True,
    type=int,
    help="Degree of the Savitsky-Golay filter polynomial.",
)
@bool_flag(
    "savgol",
    "savgol",
    True,
    "Apply Savitsky-Golay filtering to the adjusted WPS.",
)
@click.option(
    "--mean",
    "mean",
    is_flag=True,
    help="Use a mean filter instead of a median filter.",
)
@click.option(
    "--subtract-edges",
    "subtract_edges",
    is_flag=True,
    help="Subtract the mean of each interval's edges (see --edge-size) before "
    "filtering.",
)
@click.option(
    "--edge-size",
    "edge_size",
    metavar="BP",
    default=500,
    show_default=True,
    type=int,
    help="Edge width subtracted from each interval end when --subtract-edges "
    "is set.",
)
@threads_option()
@verbose_option()
def adjust_wps(**params):
    """Adjust raw WPS with a median filter and a Savitsky-Golay filter.

    \b
    INPUT        bigWig file of raw WPS results.
    REGIONS      BED file of intervals WPS was calculated over.
    CHROM_SIZES  .chrom.sizes file of chromosome names and sizes.
    """
    _dispatch(**params)


# --------------------------------------------------------------------------- #
# DELFI
# --------------------------------------------------------------------------- #
@click.command("delfi", epilog=_ex("finaletoolkit delfi sample.bam autosomes.chrom.sizes hg19.2bit bins.bed -g hg19 -o delfi.tsv"))
@input_argument()
@click.argument("chrom_sizes", metavar="CHROM_SIZES")
@click.argument("reference_file", metavar="REFERENCE")
@click.argument("bins_file", metavar="BINS")
@click.option(
    "-b",
    "--blacklist",
    "blacklist_file",
    metavar="BED",
    help="BED file of regions to ignore when calculating DELFI.",
)
@click.option(
    "-g",
    "--gap-file",
    "gap_file",
    metavar="GAPS",
    help='Telomere/centromere annotations: a BED4 file (4th column '
    '"centromere"/"telomere"/"short arm") or a genome name '
    "(hg19/b37/hg38/GRCh38).",
)
@output_option("Output file (.bed, .bed.gz, .tsv, or .csv).")
@click.option(
    "--no-gc-correct",
    "no_gc_correct",
    is_flag=True,
    default=False,
    help="Skip GC correction.",
)
@bool_flag(
    "remove-nocov",
    "remove_nocov",
    True,
    "Remove the two hg19 regions with no coverage (use --no-remove-nocov for "
    "non-hg19 reference genomes).",
)
@bool_flag(
    "merge-bins",
    "merge_bins",
    True,
    "Merge input bins to 5Mb (use --no-merge-bins to keep the input bins).",
)
@click.option(
    "--merge-size",
    "window_size",
    metavar="BP",
    default=5000000,
    show_default=True,
    type=int,
    help="Target size of merged genomic intervals.",
)
@min_mapq_option(30)
@threads_option()
@verbose_option()
def delfi(**params):
    """DELFI features: short/long fragments, DELFI ratio, and total fragments.

    \b
    INPUT        BAM/CRAM/Fragment file of fragment data.
    CHROM_SIZES  Tab-delimited chrom name and integer length (autosomes only
                 to replicate the original scripts).
    REFERENCE    A .2bit or FASTA (.fa, .fasta, .fna) reference genome file.
    BINS         BED file of bins over which to calculate DELFI.
    """
    _dispatch(**params)


# --------------------------------------------------------------------------- #
# Motifs & MDS
# --------------------------------------------------------------------------- #
@click.command("end-motifs", epilog=_ex("finaletoolkit end-motifs sample.bam hg38.2bit -k 4 -o motifs.tsv"))
@input_argument()
@click.argument("refseq_file", metavar="REFERENCE")
@kmer_option(4)
@min_length_option(50, "Minimum fragment length to include.")
@max_length_option(None, "Maximum fragment length to include.")
@strand_option("end")
@output_option("TSV of k-mer frequencies.")
@min_mapq_option(20)
@threads_option()
@verbose_option()
def end_motifs(**params):
    """Frequency of k-mer 5' end motifs.

    \b
    INPUT      BAM/CRAM/Fragment file of fragment data.
    REFERENCE  .2bit or FASTA reference genome file.
    """
    _dispatch(**params)


@click.command("interval-end-motifs", epilog=_ex("finaletoolkit interval-end-motifs sample.bam hg38.2bit intervals.bed -o interval_motifs.tsv"))
@input_argument()
@click.argument("refseq_file", metavar="REFERENCE")
@click.argument("intervals", metavar="REGIONS")
@kmer_option(4)
@min_length_option(50, "Minimum fragment length to include.")
@max_length_option(None, "Maximum fragment length to include.")
@strand_option("end")
@output_option("TSV or CSV file to write end-motif frequencies to.")
@min_mapq_option(20)
@threads_option()
@verbose_option()
def interval_end_motifs(**params):
    """Frequency of k-mer 5' end motifs in each region of a BED file.

    \b
    INPUT      BAM/CRAM/Fragment file of fragment data.
    REFERENCE  .2bit or FASTA reference genome file.
    REGIONS    BED file of intervals to retrieve end-motif frequencies over.
    """
    _dispatch(**params)


@click.command("breakpoint-motifs", epilog=_ex("finaletoolkit breakpoint-motifs sample.bam hg38.2bit -k 6 -o bp.tsv"))
@input_argument()
@click.argument("refseq_file", metavar="REFERENCE")
@kmer_option(6)
@min_length_option(50, "Minimum fragment length to include.")
@max_length_option(None, "Maximum fragment length to include.")
@strand_option("breakpoint")
@output_option("TSV of k-mer frequencies.")
@min_mapq_option(20)
@threads_option()
@verbose_option()
def breakpoint_motifs(**params):
    """Frequency of k-mer breakpoint motifs.

    \b
    INPUT      BAM/CRAM/Fragment file of fragment data.
    REFERENCE  .2bit or FASTA reference genome file.
    """
    _dispatch(**params)


@click.command("interval-breakpoint-motifs", epilog=_ex("finaletoolkit interval-breakpoint-motifs sample.bam hg38.2bit intervals.bed -o interval_bp.tsv"))
@input_argument()
@click.argument("refseq_file", metavar="REFERENCE")
@click.argument("intervals", metavar="REGIONS")
@kmer_option(6)
@min_length_option(50, "Minimum fragment length to include.")
@max_length_option(None, "Maximum fragment length to include.")
@strand_option("breakpoint")
@output_option("TSV or CSV file to write breakpoint-motif frequencies to.")
@min_mapq_option(20)
@threads_option()
@verbose_option()
def interval_breakpoint_motifs(**params):
    """Frequency of k-mer breakpoint motifs in each region of a BED file.

    \b
    INPUT      BAM/CRAM/Fragment file of fragment data.
    REFERENCE  .2bit or FASTA reference genome file.
    REGIONS    BED file of intervals to retrieve breakpoint-motif frequencies
               over.
    """
    _dispatch(**params)


@click.command("mds", epilog=_ex("finaletoolkit mds motifs.tsv"))
@click.argument("file_path", metavar="INPUT", required=False, default="-")
@click.option(
    "-s",
    "--sep",
    "sep",
    default="\t",
    show_default=False,
    help="Field separator. Default is a tab.",
)
@click.option(
    "--header",
    "header",
    default=0,
    show_default=True,
    type=int,
    help="Number of header rows to ignore.",
)
def mds(**params):
    """Motif diversity score (MDS) from k-mer frequencies (Jiang et al., 2020).

    \b
    INPUT  Tabular file with a k-mer column and a frequency column. Reads from
           stdin by default.
    """
    _dispatch(**params)


@click.command("regional-mds", epilog=_ex("finaletoolkit regional-mds interval_motifs.tsv rmds.bed"))
@click.argument("file_path", metavar="INPUT", required=False, default="-")
@click.argument("file_out", metavar="OUTPUT")
@click.option(
    "-s",
    "--sep",
    "sep",
    default="\t",
    show_default=False,
    help="Field separator. Default is a tab.",
)
@click.option(
    "--header",
    "header",
    default=0,
    show_default=True,
    type=int,
    help="Number of header rows to ignore.",
)
@click.option(
    "--miller-madow",
    "miller_madow",
    is_flag=True,
    default=False,
    help="Apply the Miller-Madow bias correction to each region's rMDS. "
    "Counteracts the downward bias of the plug-in entropy estimator in "
    "regions with few fragments.",
)
def regional_mds(**params):
    """Regional Motif Diversity Score (rMDS) for each region (Bandaru et al. 2026).

    \b
    INPUT   Tabular file of per-region k-mer frequencies. Reads from stdin by
            default.
    OUTPUT  BED/BEDGraph file of the regional MDS for each region.
    """
    _dispatch(**params)


# --------------------------------------------------------------------------- #
# Utilities
# --------------------------------------------------------------------------- #
@click.command("filter-file", epilog=_ex("finaletoolkit filter-file sample.bam -q 30 --min-length 100 -o filtered.bam"))
@click.argument("input_file", metavar="INPUT")
@reference_option()
@click.option(
    "-w",
    "--whitelist",
    "whitelist_file",
    metavar="BED",
    default=None,
    help="Only keep alignments overlapping intervals in this BED file.",
)
@click.option(
    "-b",
    "--blacklist",
    "blacklist_file",
    metavar="BED",
    default=None,
    help="Only keep alignments outside intervals in this BED file.",
)
@output_option("Output BED/BAM/CRAM file path.")
@min_mapq_option(30)
@min_length_option(None, "Minimum fragment length to include.")
@max_length_option(None, "Maximum fragment length to include.")
@intersect_policy_option()
@threads_option()
@verbose_option()
def filter_file(**params):
    """Filter a BED/BAM/CRAM file by pairing, MAPQ, flags, length, and regions.

    \b
    INPUT  BAM/CRAM/BED file to filter.
    """
    _dispatch(**params)


@click.command("agg-bw", epilog=_ex("finaletoolkit agg-bw wps.bw intervals.bed -m 120 -o agg.wig"))
@click.argument("input_file", metavar="INPUT")
@click.argument("interval_file", metavar="REGIONS")
@output_option("Wiggle file of the aggregate signal over the input intervals.")
@click.option(
    "-m",
    "--median-window-size",
    "median_window_size",
    metavar="BP",
    default=1,
    show_default=True,
    type=int,
    help="Median filter window used to aggregate scores. Set to 120 when "
    "aggregating WPS signals.",
)
@click.option(
    "--mean",
    "mean",
    is_flag=True,
    help="Use the mean instead of the sum.",
)
@verbose_option()
def agg_bw(**params):
    """Aggregate a bigWig signal over constant-length BED intervals.

    \b
    INPUT    bigWig file of signal over the input intervals.
    REGIONS  BED file of intervals the signal was calculated over.
    """
    _dispatch(**params)


@click.command("gap-bed", epilog=_ex("finaletoolkit gap-bed hg19 hg19_gaps.bed"))
@click.argument(
    "reference_genome",
    metavar="GENOME",
    type=click.Choice(["hg19", "b37", "human_g1k_v37", "hg38", "GRCh38"]),
)
@click.argument("output_file", metavar="OUTPUT")
def gap_bed(**params):
    """BED4 of centromeres, telomeres, and short-arm intervals (UCSC hg19 gaps).

    'Gap' is used liberally; for hg38/GRCh38 it may refer to regions that no
    longer have gaps in the reference.

    \b
    GENOME  Reference genome (hg19, b37, human_g1k_v37, hg38, GRCh38).
    OUTPUT  Path to write the BED file to. Pass '-' for standard output.
    """
    _dispatch(**params)


# Ordered to match the original subcommand registration order.
COMMANDS = [
    coverage,
    frag_length_bins,
    frag_length_intervals,
    cleavage_profile,
    wps,
    adjust_wps,
    delfi,
    end_motifs,
    interval_end_motifs,
    breakpoint_motifs,
    interval_breakpoint_motifs,
    mds,
    regional_mds,
    filter_file,
    agg_bw,
    gap_bed,
]


def register_commands(group: click.Group) -> None:
    """Register every FinaleToolkit subcommand on ``group``."""
    for command in COMMANDS:
        group.add_command(command)
