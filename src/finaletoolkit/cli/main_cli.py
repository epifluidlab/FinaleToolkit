#!/usr/bin/env python3
"""
Module containing command and subcommands for `finaletoolkit` CLI
standalone program.
"""

from __future__ import annotations
import argparse
from sys import stderr
import importlib

from .. import __version__


def main_cli_parser():
    """
    Returns argparse parser for CLI.
    """

    # top level parser
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

    # coverage
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
        '--output-file',
        default='-',
        help='A BED file containing coverage values over the intervals '
        'specified in interval file.')
    cli_coverage.add_argument(
        '-n',
        '--normalize',
        action='store_true',
        help="If flag set, multiplies by user inputed scale factor if"
        " given and normalizes output by total coverage. May lead to "
        "longer execution time for high-throughput data."
    )
    cli_coverage.add_argument(
        '-s',
        '--scale-factor',
        default=1.,
        type=float,
        help='Scale factor for coverage values. Default is 1.')
    cli_coverage.add_argument(
        '-min',
        '--min-length',
        default=0,
        type=int,
        help='Minimum length for a fragment to be included in coverage.'
        )
    cli_coverage.add_argument(
        '-max',
        '--max-length',
        default=None,
        type=int,
        help='Maximum length for a fragment to be included in coverage.'
        )
    cli_coverage.add_argument(
        '-p',
        '--intersect-policy',
        choices=['midpoint', 'any'],
        default='midpoint',
        type=str,
        help='Specifies what policy is used to include fragments in the'
        ' given interval. See User Guide for more information.')
    cli_coverage.add_argument(
        '-q',
        '--quality-threshold',
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
        module='finaletoolkit.frag._coverage', func='coverage')

    # frag-length-bins
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
        '-E',
        '--stop',
        help='Specify the 1-based right-most coordinate of the interval'
        ' to select fragments from. (Must also specify --contig.)',
        type=int)
    cli_frag_length_bins.add_argument(
        '-min',
        '--min-length',
        default=0,
        type=int,
        help='Minimum length for a fragment to be included in fragment length.'
        )
    cli_frag_length_bins.add_argument(
        '-max',
        '--max-length',
        default=None,
        type=int,
        help='Maximum length for a fragment to be included in fragment length.'
        )
    cli_frag_length_bins.add_argument(
        '-p',
        '--intersect-policy',
        choices=['midpoint', 'any'],
        default='midpoint',
        type=str,
        help='Specifies what policy is used to include fragments in the'
        ' given interval. See User Guide for more information.')
    cli_frag_length_bins.add_argument(
        '--bin-size',
        type=int,
        default=1,
        help='Specify the size of the bins to group fragment lengths '
        'into.')
    cli_frag_length_bins.add_argument(
        '-o',
        '--output-file',
        default='-',
        type=str,
        help='A .TSV file containing containing fragment lengths binned'
        ' according to the specified bin size.')
    cli_frag_length_bins.add_argument(
        '--histogram-path',
        default=None,
        help='Path to store histogram.',
    )
    cli_frag_length_bins.add_argument(
        '-q',
        '--quality-threshold',
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
    cli_frag_length_bins.set_defaults(module='finaletoolkit.frag._frag_length', func='frag_length_bins')

    # frag-length-intervals
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
        '-min',
        '--min-length',
        default=0,
        type=int,
        help='Minimum length for a fragment to be included in fragment length.'
        )
    cli_frag_length_intervals.add_argument(
        '-max',
        '--max-length',
        default=None,
        type=int,
        help='Maximum length for a fragment to be included in fragment length.'
        )
    cli_frag_length_intervals.add_argument(
        '-p',
        '--intersect-policy',
        choices=['midpoint', 'any'],
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
    cli_frag_length_intervals.set_defaults(module='finaletoolkit.frag._frag_length', func='frag_length_intervals')

    # cleavage-profile
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
        '-c',
        '--chrom-sizes',
        help='A .chrom.sizes file containing chromosome names and sizes.')
    cli_cleavage_profile.add_argument(
        '-o',
        '--output-file',
        default='-',
        help='A bigWig file containing the cleavage proportion results over '
        'the intervals specified in interval file.',)
    cli_cleavage_profile.add_argument(
        '-min',
        '--min-length',
        default=0,
        type=int,
        help='Minimum length for a fragment to be included.'
        )
    cli_cleavage_profile.add_argument(
        '-max',
        '--max-length',
        default=None,
        type=int,
        help='Maximum length for a fragment to be included.'
        )
    cli_cleavage_profile.add_argument(
        '-lo',
        '--fraction_low',
        type=int,
        dest='min_length',
        help="Minimum length for a fragment to be included in cleavage "
        "proportion calculation. Deprecated. Use --min-length instead.")
    cli_cleavage_profile.add_argument(
        '-hi',
        '--fraction-high',
        type=int,
        dest='max_length',
        help="Maximum length for a fragment to be included in cleavage "
        "proportion calculation. Deprecated. Use --max-length instead.")
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
        'coordinates. Default is 0.')
    cli_cleavage_profile.add_argument(
        '-r', '--right',
        default=0,
        type=int,
        help='Number of base pairs to add to stop coordinate to create '
        'interval. Useful when dealing with BED files with only CpG '
        'coordinates. Default is 0.')
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
    cli_cleavage_profile.set_defaults(module='finaletoolkit.frag._cleavage_profile', func='multi_cleavage_profile')

    # wps
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
        help='Path to a BED file containing sites to calculate WPS '
        'over. The intervals in this BED file should be sorted, first '
        'by `contig` then `start`.')
    cli_wps.add_argument(
        '-c',
        '--chrom-sizes',
        help='A .chrom.sizes file containing chromosome names and sizes.')
    cli_wps.add_argument(
        '-o',
        '--output-file',
        default='-',
        help='A bigWig file containing the WPS results over the '
        'intervals specified in interval file.')
    cli_wps.add_argument(
        '-i',
        '--interval-size',
        default=5000,
        type=int,
        help='Size in bp of the intervals to calculate WPS over. These'
        'new intervals are centered over those specified in the site_bed.'
        'Default is 5000')
    cli_wps.add_argument(
        '-W',
        '--window-size',
        default=120,
        type=int,
        help='Size of the sliding window used to calculate WPS scores.'
        ' Default is 120')
    cli_wps.add_argument(
        '-min',
        '--min-length',
        default=0,
        type=int,
        help='Minimum length for a fragment to be included.'
        )
    cli_wps.add_argument(
        '-max',
        '--max-length',
        default=None,
        type=int,
        help='Maximum length for a fragment to be included.'
        )
    cli_wps.add_argument(
        '-lo',
        '--fraction_low',
        type=int,
        dest='min_length',
        help='Minimum length for a fragment to be included in WPS '
        'calculation. Deprecated. Use --min-length instead.')
    cli_wps.add_argument(
        '-hi',
        '--fraction_high',
        type=int,
        dest='max_length',
        help='Maximum length for a fragment to be included in WPS '
        'calculation. Deprecated. Use --max-length instead.')
    cli_wps.add_argument(
        '-q',
        '--quality-threshold',
        default=30,
        type=int,
        help="Minimum mapping quality threshold. Default is 30")
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
    cli_wps.set_defaults(module='finaletoolkit.frag', func='multi_wps')

    # adjust wps
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
        'chrom_sizes',
        help='A .chrom.sizes file containing chromosome names and sizes.')
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
        help='Size of the median filter or mean filter window used to adjust WPS scores.')
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
        '-S',
        '--exclude-savgol',
        dest='savgol',
        action='store_false',
        help='Do not perform Savitsky-Golay filtering'
        'scores.')
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
        '--edge-size',
        default=500,
        help='size of the edge subtracted from ends of window when '
        '--subtract-edges is set. Default is 500.')
    cli_adjust_wps.add_argument(
        '-v',
        '--verbose',
        action='count',
        help='Enable verbose mode to display detailed processing information.')
    cli_adjust_wps.set_defaults(module='finaletoolkit.frag._adjust_wps', func='adjust_wps')

    # delfi
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
        'chrom_sizes',
        help="Tab-delimited file containing (1) chrom name and (2) integer "
        "length of chromosome in base pairs. Should contain only autosomes if"
        "You want to replicate the original scripts.")
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
        '-s',
        '--window-size',
        default=5000000,    # TODO: consider accepting 5Mb or similar format.
        help='Specify size of large genomic intervals to merge smaller 100kb '
        'intervals (or whatever the user specified in bins_file) into. Default'
        'is 5000000'
    )
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
    cli_delfi.set_defaults(module='finaletoolkit.frag._delfi', func='delfi')

    # delfi-gc-correct
    cli_delfi_gc = subparsers.add_parser(
        'delfi-gc-correct',
        prog='finaletoolkit-delfi-gc-correct',
        description='Performs gc-correction on raw delfi data. This '
        'command is deprecated and will be removed in a future version '
        'of FinaleToolkit. The delfi command has gc correction on by '
        'default.')
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
    cli_delfi_gc.set_defaults(module='finaletoolkit.frag._delfi_gc_correct', func='cli_delfi_gc_correct')

    # end-motifs
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
        '-min',
        '--min-length',
        default=0,
        type=int,
        help='Minimum length for a fragment to be included.'
        )
    cli_motifs.add_argument(
        '-max',
        '--max-length',
        default=None,
        type=int,
        help='Maximum length for a fragment to be included.'
        )
    cli_motifs.add_argument(
        '-B',
        '--no-both-strands',
        action="store_false",
        dest="both_strands",
        help="Set flag to only consider one strand for end-motifs."
    )
    cli_motifs.add_argument(
        '-n',
        '--negative-strand',
        action="store_true",
        help="Set flag in conjunction with -B to only consider 5' end motifs "
        "on the negative strand."
    )
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
    cli_motifs.set_defaults(module='finaletoolkit.frag._end_motifs', func='end_motifs')

    # interval-end-motifs
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
        '-min',
        '--min-length',
        default=0,
        type=int,
        help='Minimum length for a fragment to be included.'
        )
    cli_interval_motifs.add_argument(
        '-max',
        '--max-length',
        default=None,
        type=int,
        help='Maximum length for a fragment to be included.'
        )
    cli_interval_motifs.add_argument(
        '-lo',
        '--fraction-low',
        type=int,
        dest='min_length',
        help='Deprecated alias for --min-length')
    cli_interval_motifs.add_argument(
        '-hi',
        '--fraction-high',
        type=int,
        dest='max_length',
        help='Deprecated alias for --max-length')
    cli_interval_motifs.add_argument(
        '-B',
        '--single-strand',
        action="store_false",
        dest="both_strands",
        help="Set flag to only consider one strand for end-motifs. By default,"
        " the positive strand is calculated, but with the -n flag, the "
        "5' end motifs of the negative strand are considered instead."
    )
    cli_interval_motifs.add_argument(
        '-n',
        '--negative-strand',
        action="store_true",
        help="Set flag in conjunction with -B to only consider 5' end motifs "
        "on the negative strand."
    )
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
    cli_interval_motifs.set_defaults(module='finaletoolkit.frag._end_motifs', func='interval_end_motifs')

    # mds
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
    cli_mds.set_defaults(module='finaletoolkit.frag._end_motifs', func='_cli_mds')

    # interval-mds
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
    cli_interval_mds.add_argument(
        '--header',
        default=0,
        type=int,
        help='Number of header rows to ignore. Default is 0')
    cli_interval_mds.set_defaults(module='finaletoolkit.frag._end_motifs', func='_cli_interval_mds')

    # filter-bam
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
        '--quality-threshold',
        type=int,
        default=30,
        help='Minimum mapping quality threshold.')
    cli_filter_bam.add_argument(
        '-min',
        '--min-length',
        default=None,
        type=int,
        help='Minimum length for a fragment to be included.'
        )
    cli_filter_bam.add_argument(
        '-max',
        '--max-length',
        default=None,
        type=int,
        help='Maximum length for a fragment to be included.'
        )
    cli_filter_bam.add_argument(
        '-lo',
        '--fraction-low',
        type=int,
        dest='min_length',
        help='Deprecated alias for --min-length')
    cli_filter_bam.add_argument(
        '-hi',
        '--fraction-high',
        type=int,
        dest='max_length',
        help='Deprecated alias for --max-length')
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
    cli_filter_bam.set_defaults(module='finaletoolkit.utils', func='filter_bam')

    # agg-bw
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
        default=1,
        type=int,
        help='Size of the median filter window used to aggregate '
        'scores. Set to 120 if aggregating WPS signals.')
    cli_agg_bw.add_argument(
        '-a',
        '--mean',
        action='store_true',
        help='use mean instead'
    )
    cli_agg_bw.add_argument(
        '-v',
        '--verbose',
        action='count',
        help='Enable verbose mode to display detailed processing '
        'information.')
    cli_agg_bw.set_defaults(module='finaletoolkit.utils._agg_bw', func='agg_bw')

    # gap-bed
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
        choices=['hg19', 'b37', 'human_g1k_v37', 'hg38', 'GRCh38'],
        help='Reference genome to provide gaps for.')
    cli_gap_bed.add_argument(
        'output_file',
        help='Path to write BED file to. If "-" used, writes to stdout.')
    cli_gap_bed.set_defaults(module='finaletoolkit.genome.gaps', func='_cli_gap_bed')

    return parser


def main_cli():
    """
    Function called when using `finaletoolkit` CLI standalone program.
    """
    parser = main_cli_parser()

    args = parser.parse_args()
    if hasattr(args, "func"):
        try:
            # lazy loading function
            funcargs = vars(args)
            func_module = funcargs.pop('module')
            func_name = funcargs.pop('func')
            
            module = importlib.import_module(func_module)
            function = getattr(module, func_name)

            function(**funcargs)
        except AttributeError as e:
            stderr.write(f"FinaleToolkit recieved AttributeError: {e}\n")
            stderr.write("Please see usage instructions below.\n")
            parser.print_help()
    else:
        parser.print_help()


if __name__ == '__main__':
    main_cli()
