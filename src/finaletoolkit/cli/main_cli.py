#!/usr/bin/env python3
"""
Module containing command and subcommands for `finaletoolkit` CLI
standalone program.
"""

from __future__ import annotations
import argparse

from finaletoolkit import __version__
from finaletoolkit.frag.frag_length import (
    frag_length_bins, frag_length_intervals
)
from finaletoolkit.utils.agg_bw import agg_bw
from finaletoolkit.utils.filter_bam import filter_bam
from finaletoolkit.frag.coverage import coverage
from finaletoolkit.frag.multi_wps import multi_wps
from finaletoolkit.frag.delfi import delfi
from finaletoolkit.frag.adjust_wps import adjust_wps
from finaletoolkit.frag.delfi_gc_correct import cli_delfi_gc_correct
from finaletoolkit.frag.end_motifs import (
    end_motifs, _cli_mds, _cli_interval_mds, interval_end_motifs)
from finaletoolkit.frag.cleavage_profile import _cli_cleavage_profile
from finaletoolkit.genome.gaps import _cli_gap_bed


def main_cli_parser():
    """
    Returns argparse parser for CLI.
    """

    parser = argparse.ArgumentParser(
        description='FinaleToolkit is a package and standalone program '
        'to extract fragmentation features of cell-free DNA from '
        'paired-end sequencing data.',
        epilog='')
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=f'FinaleToolkit {__version__}')
    subparsers = parser.add_subparsers()

    cli_coverage = subparsers.add_parser(
        'coverage',
        description='Calculates fragmentation coverage over intervals '
        'defined in a BED file based on alignment data from a '
        'BAM/SAM/CRAM/Fragment file.')
    cli_coverage.add_argument(
        'input_file',
        help='Path to a BAM/SAM/CRAM/Fragment file containing fragment '
        'data.')
    cli_coverage.add_argument(
        'interval_file',
        help='Path to a BED file containing intervals to calculate '
        'coverage over.')
    cli_coverage.add_argument(
        '-o',
        '--output_file',
        default='-',
        help='A BED file containing coverage values over the intervals '
        'specified in interval file.')
    cli_coverage.add_argument(
        '-n',
        '--normalize',
        action='store_true',
        help="If flag set, ignores any user inputed scale factor and "
        "normalizes output by total coverage."
    )
    cli_coverage.add_argument(
        '-s',
        '--scale-factor',
        default=1e6,
        type=float,
        help='Scale factor for coverage values.')
    cli_coverage.add_argument(
        '-p',
        '--intersect_policy',
        choices=['midpoint',
        'any'],
        default='midpoint',
        type=str,
        help='Specifies what policy is used to include fragments in the'
        ' given interval. See User Guide for more information.')
    cli_coverage.add_argument(
        '-q',
        '--quality_threshold',
        default=30,
        type=int,
        help='Minimum mapping quality threshold.')
    cli_coverage.add_argument(
        '-w',
        '--workers',
        default=1,
        type=int,
        help='Number of worker processes.')
    cli_coverage.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        help='Enable verbose mode to display detailed processing information.')
    cli_coverage.set_defaults(
        func=coverage)

    cli_frag_length_bins = subparsers.add_parser(
        'frag-length-bins',
        prog='finaletoolkit-frag-length-bins',
        description='Retrieves fragment lengths grouped in bins given a'
        ' BAM/SAM/CRAM/Fragment file.')
    cli_frag_length_bins.add_argument(
        'input_file',
        help='Path to a BAM/SAM/CRAM/Fragment file containing fragment data.')
    cli_frag_length_bins.add_argument(
        '-c',
        '--contig',
        type=str,
        help='Specify the contig or chromosome to select fragments from'
        '. (Required if using --start or --stop.)')
    cli_frag_length_bins.add_argument(
        '-S',
        '--start',
        type=int,
        help='Specify the 0-based left-most coordinate of the interval '
        'to select fragments from. (Must also specify --contig.)')
    cli_frag_length_bins.add_argument(
        '-p',
        '--intersect_policy',
        choices=['midpoint',
        'any'],
        default='midpoint',
        type=str,
        help='Specifies what policy is used to include fragments in the'
        ' given interval. See User Guide for more information.')
    cli_frag_length_bins.add_argument(
        '-E',
        '--stop',
        help='Specify the 1-based right-most coordinate of the interval'
        ' to select fragments from. (Must also specify --contig.)',
        type=int)
    cli_frag_length_bins.add_argument(
        '--bin-size',
        type=int,
        help='Specify the size of the bins to group fragment lengths '
        'into.')
    cli_frag_length_bins.add_argument(
        '-o',
        '--output_file',
        default='-',
        type=str,
        help='A .TSV file containing containing fragment lengths binned'
        ' according to the specified bin size.')
    cli_frag_length_bins.add_argument(
        '--contig-by-contig',
        action='store_true',
        help='Placeholder, not implemented.')
    cli_frag_length_bins.add_argument(
        '--histogram',
        action='store_true',
        help='Enable histogram mode to display histogram in terminal.')
    cli_frag_length_bins.add_argument(
        '-q',
        '--quality_threshold',
        default=30,
        type=int,
        help="Minimum mapping quality threshold.")
    cli_frag_length_bins.add_argument(
        '-v',
        '--verbose',
        action='count',
        default=0,
        help='Enable verbose mode to display detailed processing '
        'information.')
    cli_frag_length_bins.set_defaults(func=frag_length_bins)

    cli_frag_length_intervals = subparsers.add_parser(
        'frag-length-intervals',
        description='Retrieves fragment length summary statistics over '
        'intervals defined in a BED file based on alignment data from a'
        ' BAM/SAM/CRAM/Fragment file.')
    cli_frag_length_intervals.add_argument(
        'input_file',
        help='Path to a BAM/SAM/CRAM/Fragment file containing fragment '
        'data.')
    cli_frag_length_intervals.add_argument(
        'interval_file',
        help='Path to a BED file containing intervals to retrieve '
        'fragment length summary statistics over.')
    cli_frag_length_intervals.add_argument(
        '-p',
        '--intersect_policy',
        choices=['midpoint',
        'any'],
        default='midpoint',
        type=str,
        help='Specifies what policy is used to include fragments in the'
        ' given interval. See User Guide for more information.')
    cli_frag_length_intervals.add_argument(
        '-o',
        '--output-file',
        default='-',
        help='A BED file containing fragment length summary statistics '
        '(mean, median, st. dev, min, max) over the intervals specified'
        ' in the interval file.')
    cli_frag_length_intervals.add_argument(
        '-q',
        '--quality-threshold',
        default=30,
        type=int,
        help='Minimum mapping quality threshold.')
    cli_frag_length_intervals.add_argument(
        '-w',
        '--workers',
        default=1,
        type=int,
        help='Number of worker processes.')
    cli_frag_length_intervals.add_argument(
        '-v',
        '--verbose',
        default=0,
        action='count',
        help='Enable verbose mode to display detailed processing '
        'information.')
    cli_frag_length_intervals.set_defaults(func=frag_length_intervals)

    cli_cleavage_profile = subparsers.add_parser(
        'cleavage-profile',
        prog='finaletoolkit-cleavage-profile',
        description='Calculates cleavage proportion over intervals defined in '
        'a BED file based on alignment data from a BAM/SAM/CRAM/Fragment '
        'file.')
    cli_cleavage_profile.add_argument(
        'input_file',
        help='Path to a BAM/SAM/CRAM/Fragment file containing fragment data.')
    cli_cleavage_profile.add_argument(
        'interval_file',
        help='Path to a BED file containing intervals to calculates cleavage '
        'proportion over.')
    cli_cleavage_profile.add_argument(
        '-o',
        '--output_file',
        default='-',
        help='A bigWig file containing the cleavage proportion results over '
        'the intervals specified in interval file.',)
    cli_cleavage_profile.add_argument(
        '-lo',
        '--fraction_low',
        default=120,
        type=int,
        help="Minimum length for a fragment to be included in cleavage "
        "proportion calculation.")
    cli_cleavage_profile.add_argument(
        '-hi',
        '--fraction_high',
        default=180,
        type=int,
        help="Maximum length for a fragment to be included in cleavage "
        "proportion calculation.")
    cli_cleavage_profile.add_argument(
        '-q',
        '--quality-threshold',
        default=20,
        type=int,
        help='Minimum mapping quality threshold.')
    cli_cleavage_profile.add_argument(
        '-l',
        '--left',
        default=0,
        type=int,
        help='Number of base pairs to subtract from start coordinate to create'
        ' interval. Useful when dealing with BED files with only CpG '
        'coordinates.')
    cli_cleavage_profile.add_argument(
        '-r','--right',
        default=0,
        type=int,
        help='Number of base pairs to add to stop coordinate to create '
        'interval. Useful when dealing with BED files with only CpG '
        'coordinates.')
    cli_cleavage_profile.add_argument(
        '-w',
        '--workers',
        default=1,
        type=int,
        help='Number of worker processes.')
    cli_cleavage_profile.add_argument(
        '-v',
        '--verbose',
        default=0,
        action='count',
        help='Enable verbose mode to display detailed processing information.')
    cli_cleavage_profile.set_defaults(func=_cli_cleavage_profile)

    cli_wps = subparsers.add_parser(
        'wps',
        prog='finaletoolkit-wps',
        description='Calculates Windowed Protection Score (WPS) over '
        'intervals defined in a BED file based on alignment data from a'
        ' BAM/SAM/CRAM/Fragment file.')
    cli_wps.add_argument(
        'input_file',
        help='Path to a BAM/SAM/CRAM/Fragment file containing fragment '
        'data.')
    cli_wps.add_argument(
        'site_bed',
        help='Path to a BED file containing intervals to calculate WPS '
        'over. The intervals in this BED file should be sorted, first '
        'by `contig` then `start`.')
    cli_wps.add_argument(
        '-o',
        '--output_file',
        default='-',
        help='A bigWig file containing the WPS results over the '
        'intervals specified in interval file.')
    cli_wps.add_argument(
        '-i',
        '--interval_size',
        default=5000,
        type=int,
        help='Size in bp of each interval in the interval file.')
    cli_wps.add_argument(
        '-W',
        '--window_size',
        default=120,
        type=int,
        help='Size of the sliding window used to calculate WPS scores.')
    cli_wps.add_argument(
        '-lo',
        '--fraction_low',
        default=120,
        type=int,
        help='Minimum length for a fragment to be included in WPS '
        'calculation.')
    cli_wps.add_argument(
        '-hi',
        '--fraction_high',
        default=180,
        type=int,
         help='Maximum length for a fragment to be included in WPS '
         'calculation.')
    cli_wps.add_argument(
        '-q',
        '--quality_threshold',
        default=30,
        type=int,
        help="Minimum mapping quality threshold.")
    cli_wps.add_argument(
        '-w',
        '--workers',
        default=1,
        type=int,
        help='Number of worker processes.')
    cli_wps.add_argument(
        '-v',
        '--verbose',
        action='count',
        default=0,
        help='Enable verbose mode to display detailed processing information.')
    cli_wps.set_defaults(func=multi_wps)

    cli_adjust_wps = subparsers.add_parser(
        'adjust-wps',
        prog='finaletoolkit-adjust-wps',
        description='Adjusts raw Windowed Protection Score (WPS) by applying a'
        ' median filter and Savitsky-Golay filter.')
    cli_adjust_wps.add_argument(
        'input_file',
        help='A bigWig file containing the WPS results over the intervals '
        'specified in interval file.')
    cli_adjust_wps.add_argument(
        'interval_file',
        help='Path to a BED file containing intervals to WPS was calculated '
        'over.')
    cli_adjust_wps.add_argument(
        'genome_file',
        help='A .chrom.sizes file containing chromosome sizes.')
    cli_adjust_wps.add_argument(
        '-o',
        '--output-file',
        default='-',
        help='A bigWig file containing the adjusted WPS results over the '
        'intervals specified in interval file.')
    cli_adjust_wps.add_argument(
        '-i',
        '--interval_size',
        default=5000,
        type=int,
        help='Size in bp of each interval in the interval file.')
    cli_adjust_wps.add_argument(
        '-m',
        '--median-window-size',
        default=1000,
        type=int,
        help='Size of the median filter window used to adjust WPS scores.')
    cli_adjust_wps.add_argument(
        '-s',
        '--savgol-window-size',
        default=21,
        type=int,
        help='Size of the Savitsky-Golay filter window used to adjust WPS '
        'scores.')
    cli_adjust_wps.add_argument(
        '-p',
        '--savgol-poly-deg',
        default=2,
        type=int,
        help='Degree polynomial for Savitsky-Golay filter.')
    cli_adjust_wps.add_argument(
        '-w',
        '--workers',
        default=1,
        type=int,
        help='Number of worker processes.')
    cli_adjust_wps.add_argument(
        '--mean',
        action='store_true',
        help='A mean filter is used instead of median.')
    cli_adjust_wps.add_argument(
        '--subtract-edges',
        action='store_true',
        help='Take the median of the first and last 500 bases in a window and '
        'subtract from the whole interval.')
    cli_adjust_wps.add_argument(
        '-v',
        '--verbose',
        action='count',
        help='Enable verbose mode to display detailed processing information.')
    cli_adjust_wps.set_defaults(func=adjust_wps)

    cli_delfi = subparsers.add_parser(
        'delfi',
        prog='finaletoolkit-delfi',
        description='Calculates DELFI features over genome, returning '
        'information about (GC-corrected) short fragments, long '
        'fragments, DELFI ratio, and total fragments.'
        )
    cli_delfi.add_argument(
        'input_file',
        help="Path to a BAM/SAM/CRAM/Fragment file containing fragment data.")
    cli_delfi.add_argument(
        'autosomes',
        help="Tab-delimited file containing (1) autosome name and (2) integer "
        "length of chromosome in base pairs.")
    cli_delfi.add_argument(
        'reference_file',
        help="The .2bit file for the associate reference genome sequence used "
        "during alignment.")
    cli_delfi.add_argument(
        'bins_file',
        help="A BED file containing bins over which to calculate DELFI. To "
        "replicate Cristiano et al.'s methodology, use 100kb bins over human "
        "autosomes.")
    cli_delfi.add_argument(
        '-b',
        '--blacklist-file',
        help="BED file containing regions to ignore when calculating DELFI.")
    cli_delfi.add_argument(
        '-g',
        '--gap-file',
        help='BED4 format file containing columns for "chrom", "start",'
        '"stop", and "type". The "type" column should denote whether '
        'the entry corresponds to a "centromere", "telomere", or "short'
        ' arm", and entries not falling into these categories are '
        'ignored. This information corresponds to the "gap" track for '
        'hg19 in the UCSC Genome Browser.')
    cli_delfi.add_argument(
        '-o',
        '--output-file',
        default='-',
        help='BED, bed.gz, TSV, or CSV file to write DELFI data to. If '
        '"-", writes to stdout.')
    cli_delfi.add_argument(
        '-G',
        '--no-gc-correct',
        action='store_false',
        dest="gc_correct",
        help="Skip GC correction.")
    cli_delfi.add_argument(
        '-R',
        '--keep-nocov',
        action='store_false',
        dest="remove_nocov",
        help="Skip removal two regions in hg19 with no coverage. Use this flag"
        " when not using hg19 human reference genome.")
    cli_delfi.add_argument(
        '-M',
        '--no-merge-bins',
        action='store_false',
        dest="merge_bins",
        help="Keep 100kb bins and do not merge to 5Mb size.")
    cli_delfi.add_argument(
        '-q',
        '--quality-threshold',
        default=30,
        type=int,
        help="Minimum mapping quality threshold.")
    cli_delfi.add_argument(
        '-w',
        '--workers',
        default=1,
        type=int,
        help="Number of worker processes.")
    cli_delfi.add_argument(
        '-v',
        '--verbose',
        action='count',
        default=0,
        help='Enable verbose mode to display detailed processing information.')
    cli_delfi.set_defaults(func=delfi)

    cli_delfi_gc = subparsers.add_parser(
        'delfi-gc-correct',
        prog='finaletoolkit-delfi-gc-correct',
        description='Performs gc-correction on raw delfi data.')
    cli_delfi_gc.add_argument(
        'input_file',
        help='BED file containing raw DELFI data. Raw DELFI data should'
        ' only have columns for "contig", "start", "stop", "arm", '
        '"short", "long", "gc", "num_frags", "ratio".')
    cli_delfi_gc.add_argument(
        '-o',
        '--output-file',
        default='-',
        help='BED to print GC-corrected DELFI fractions. If "-", will '
        'write to stdout.')
    cli_delfi_gc.add_argument(
        '--header-lines',
        default=1,
        type=int,
        help='Number of header lines in BED.')
    cli_delfi_gc.add_argument(
        '-v',
        '--verbose',
        action='count',
        help='Enable verbose mode to display detailed processing information.')
    cli_delfi_gc.set_defaults(func=cli_delfi_gc_correct)

    cli_motifs = subparsers.add_parser(
        'end-motifs',
        prog='finaletoolkit-end-motifs',
        description="Measures frequency of k-mer 5' end motifs.")
    cli_motifs.add_argument(
        'input_file',
        help='Path to a BAM/SAM/CRAM/Fragment file containing fragment data.')
    cli_motifs.add_argument(
        'refseq_file',
        help='The .2bit file for the associate reference genome sequence used '
        'during alignment.')
    cli_motifs.add_argument(
        '-k',
        default=4,
        type=int,
        help='Length of k-mer.')
    cli_motifs.add_argument(
        '-o',
        '--output-file',
        default='-',
        help='TSV to print k-mer frequencies. If "-", will write to stdout.')
    cli_motifs.add_argument(
        '-q',
        '--quality-threshold',
        default=20,
        type=int,
        help='Minimum mapping quality threshold.')
    cli_motifs.add_argument(
        '-w',
        '--workers',
        default=1,
        type=int,
        help='Number of worker processes.')
    cli_motifs.add_argument(
        '-v',
        '--verbose',
        default=0,
        action='count',
        help='Enable verbose mode to display detailed processing information.')
    cli_motifs.set_defaults(func=end_motifs)


    cli_interval_motifs = subparsers.add_parser(
        'interval-end-motifs',
        prog='finaletoolkit-interval-end-motifs',
        description="Measures frequency of k-mer 5' end motifs in each region "
        "specified in a BED file and writes data into a table.")
    cli_interval_motifs.add_argument(
        'input_file',
        help='Path to a BAM/SAM/CRAM/Fragment file containing fragment data.')
    cli_interval_motifs.add_argument(
        'refseq_file',
        help='The .2bit file for the associate reference genome sequence used '
        'during alignment.')
    cli_interval_motifs.add_argument(
        'intervals',
        help='Path to a BED file containing intervals to retrieve end motif '
        'frequencies over.')
    cli_interval_motifs.add_argument(
        '-k',
        default=4,
        type=int,
        help='Length of k-mer.')
    cli_interval_motifs.add_argument(
        '-lo',
        '--fraction-low',
        default=10,
        type=int,
        help='Minimum length for a fragment to be included in end motif '
        'frequency.')
    cli_interval_motifs.add_argument(
        '-hi',
        '--fraction-high',
        default=600,
        type=int,
        help='Maximum length for a fragment to be included in end motif '
        'frequency.')
    cli_interval_motifs.add_argument(
        '-o',
        '--output-file',
        default='-',
        help="Path to TSV or CSV file to write end motif frequencies to.")
    cli_interval_motifs.add_argument(
        '-q',
        '--quality-threshold',
        default=20,
        type=int,
        help='Minimum mapping quality threshold.')
    cli_interval_motifs.add_argument(
        '-w',
        '--workers',
        default=1,
        type=int,
        help='Number of worker processes.')
    cli_interval_motifs.add_argument(
        '-v',
        '--verbose',
        default=0,
        action='count',
        help='Enable verbose mode to display detailed processing information.')
    cli_interval_motifs.set_defaults(func=interval_end_motifs)

    cli_mds = subparsers.add_parser(
        'mds',
        prog='finaletoolkit-mds',
        description='Reads k-mer frequencies from a file and calculates a '
        'motif diversity score (MDS) using normalized Shannon entropy as '
        'described by Jiang et al (2020).')
    cli_mds.add_argument(
        'file_path',
        nargs='?',
        default='-',
        help='Tab-delimited or similar file containing one column for all '
        'k-mers a one column for frequency. Reads from stdin by default.')
    cli_mds.add_argument(
        '-s',
        '--sep',
        default='\t',
        help='Separator used in tabular file.')
    cli_mds.add_argument(
        '--header',
        default=0,
        type=int,
        help='Number of header rows to ignore. Default is 0')
    cli_mds.set_defaults(func=_cli_mds)

    cli_interval_mds = subparsers.add_parser(
        'interval-mds',
        prog='finaletoolkit-interval-mds',
        description='Reads k-mer frequencies from a file and calculates a '
        'motif diversity score (MDS) for each interval using normalized '
        'Shannon entropy as described by Jiang et al (2020).')
    cli_interval_mds.add_argument(
        'file_path',
        nargs='?',
        default='-',
        help='Tab-delimited or similar file containing one column for all '
        'k-mers a one column for frequency. Reads from stdin by default.')
    cli_interval_mds.add_argument(
        '-s',
        '--sep',
        default='\t',
        help='Separator used in tabular file.')
    cli_interval_mds.add_argument(
        'file_out',
        default='-',
        help='Path to the output BED/BEDGraph file containing MDS for each '
        'interval.')
    cli_interval_mds.set_defaults(func=_cli_interval_mds)

    cli_filter_bam = subparsers.add_parser(
        'filter-bam',
        prog='finaletoolkit-filter-bam',
        description='Filters a BAM file so that all reads are in mapped pairs,'
        ' exceed a certain MAPQ, are not flagged for quality, are read1, are '
        'not secondary or supplementary alignments, and are on the same '
        'reference sequence as the mate.')
    cli_filter_bam.add_argument(
        'input_file',
        help='Path to BAM file.')
    cli_filter_bam.add_argument(
        '-r',
        '--region-file',
        default=None,
        help='Only output alignments overlapping the intervals in this BED '
        'file will be included.')
    cli_filter_bam.add_argument(
        '-o',
        '--output-file',
        default='-',
        help='Output BAM file path.')
    cli_filter_bam.add_argument(
        '-q',
        '--quality_threshold',
        type=int,
        default=30,
        help='Minimum mapping quality threshold.')
    cli_filter_bam.add_argument(
        '-hi',
        '--fraction-high',
        type=int,
        default=None,
        help='Maximum length for a fragment to be included in output BAM.')
    cli_filter_bam.add_argument(
        '-lo',
        '--fraction-low',
        type=int,
        default=None,
        help='Minimum length for a fragment to be included in output BAM.')
    cli_filter_bam.add_argument(
        '-w',
        '--workers',
        type=int,
        default=1,
        help='Number of worker processes.')
    cli_filter_bam.add_argument(
        '-v',
        '--verbose',
        action='count',
        help='Enable verbose mode to display detailed processing information.')
    cli_filter_bam.set_defaults(func=filter_bam)

    cli_agg_bw = subparsers.add_parser(
        'agg-bw',
        prog='finaletoolkit-agg-wps',
        description='Aggregates a bigWig signal over constant-length intervals'
        ' defined in a BED file.')
    cli_agg_bw.add_argument(
        'input_file',
        help=' A bigWig file containing signals over the intervals specified '
        'in interval file.')
    cli_agg_bw.add_argument(
        'interval_file',
        help='Path to a BED file containing intervals over which signals were '
        'calculated over.')
    cli_agg_bw.add_argument(
        '-o',
        '--output-file',
        default='-',
        help='A wiggle file containing the aggregate signal over the intervals'
        ' specified in interval file.')
    cli_agg_bw.add_argument(
        '-m',
        '--median-window-size',
        default=0,
        type=int,
        help='Size of the median filter window used to adjust WPS '
        'scores. Only modify if aggregating WPS signals.')
    cli_agg_bw.add_argument(
        '-v',
        '--verbose',
        action='count',
        help='Enable verbose mode to display detailed processing '
        'information.')
    cli_agg_bw.set_defaults(func=agg_bw)

    cli_gap_bed = subparsers.add_parser(
        'gap-bed',
        prog='finaletoolkit-gap-bed',
        description='Creates a BED4 file containing centromeres, telomeres, '
        'and short-arm intervals, similar to the gaps annotation track for '
        'hg19 found on the UCSC Genome Browser (Kent et al 2002). Currently '
        'only supports hg19, b37, human_g1k_v37, hg38, and GRCh38',
        epilog='Gap is used liberally in this command, and in the case '
        'hg38/GRCh38, may refer to regions where there no longer are gaps in '
        'the reference sequence.')
    cli_gap_bed.add_argument(
        'reference_genome',
        choices=['hg19',
        'b37','human_g1k_v37',
        'hg38',
        'GRCh38'],
        help='Reference genome to provide gaps for.')
    cli_gap_bed.add_argument(
        'output_file',
        help='Path to write BED file to. If "-" used, writes to stdout.')
    cli_gap_bed.set_defaults(func=_cli_gap_bed)

    return parser


def main_cli():
    """
    Function called when using `finaletoolkit` CLI standalone program.
    """
    parser = main_cli_parser()

    args = parser.parse_args()
    try:
        function = args.func
        funcargs = vars(args)
        funcargs.pop('func')

        function(**funcargs)
    except AttributeError:
        parser.print_help()
    
if __name__ == '__main__':
    main_cli()
