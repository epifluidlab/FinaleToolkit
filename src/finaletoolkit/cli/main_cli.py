#!/usr/bin/env python3

from __future__ import annotations
import argparse

from finaletoolkit.frag.frag_length import (
    _cli_frag_length, frag_length_bins, frag_length_intervals
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
    parser = argparse.ArgumentParser(description='FinaleToolkit is a package and standalone program to extract fragmentation features of cell-free DNA from paired-end sequencing data.', epilog='')
    subparsers = parser.add_subparsers()

    parser_command1 = subparsers.add_parser('coverage', description='Calculates fragmentation coverage over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.')
    parser_command1.add_argument('input_file', help='Path to a BAM/SAM/CRAM/Fragment file containing fragment data.')
    parser_command1.add_argument('interval_file', help='Path to a BED file containing intervals to calculate coverage over.')
    parser_command1.add_argument('-o', '--output_file', default='-', help='A BED file containing coverage values over the intervals specified in interval file.')
    parser_command1.add_argument('-s', '--scale-factor', default=1e6, type=float, help='Scale factor for coverage values.')
    parser_command1.add_argument('-q', '--quality_threshold', default=30, type=int, help='Minimum mapping quality threshold.')
    parser_command1.add_argument('-w', '--workers', default=1, type=int, help='Number of worker processes.')
    parser_command1.add_argument('-v', '--verbose', action='store_true', help='Enable verbose mode to display detailed processing information.')
    parser_command1.set_defaults(func=coverage)

    parser_command3 = subparsers.add_parser('frag-length-bins', prog='finaletoolkit-frag-length-bins', description='Retrieves fragment lengths grouped in bins given a BAM/SAM/CRAM/Fragment file.')
    parser_command3.add_argument('input_file', help='Path to a BAM/SAM/CRAM/Fragment file containing fragment data.')
    parser_command3.add_argument('-c', '--contig', type=str, help='Specify the contig or chromosome to select fragments from. (Required if using --start or --stop.)')
    parser_command3.add_argument('-S', '--start', type=int, help='Specify the 0-based left-most coordinate of the interval to select fragments from. (Must also specify --contig.)')
    parser_command3.add_argument('-p', '--intersect_policy', choices=['midpoint', 'any'], default='midpoint', type=str, help='Specifies what policy is used to include fragments in the given interval. See User Guide for more information.')
    parser_command3.add_argument('-E', '--stop', help='Specify the 1-based right-most coordinate of the interval to select fragments from. (Must also specify --contig.)', type=int)
    parser_command3.add_argument('--bin-size', type=int, help='Specify the size of the bins to group fragment lengths into.')
    parser_command3.add_argument('-o', '--output_file', default='-', type=str, help='A .TSV file containing containing fragment lengths binned according to the specified bin size.')
    parser_command3.add_argument('--contig-by-contig', action='store_true', help='Placeholder, not implemented.')
    parser_command3.add_argument('--histogram', action='store_true', help='Enable histogram mode to display histogram in terminal.')
    parser_command3.add_argument('-q', '--quality_threshold', default=30, type=int, help="Minimum mapping quality threshold.")
    parser_command3.add_argument('-v', '--verbose', action='count', default=0, help='Enable verbose mode to display detailed processing information.')
    parser_command3.set_defaults(func=frag_length_bins)

    parser_command3_1 = subparsers.add_parser('frag-length-intervals', description='Retrieves fragment length summary statistics over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.')
    parser_command3_1.add_argument('input_file', help='Path to a BAM/SAM/CRAM/Fragment file containing fragment data.')
    parser_command3_1.add_argument('interval_file', help='Path to a BED file containing intervals to retrieve fragment length summary statistics over.')
    parser_command3_1.add_argument('-p', '--intersect_policy', choices=['midpoint', 'any'], default='midpoint', type=str, help='Specifies what policy is used to include fragments in the given interval. See User Guide for more information.')
    parser_command3_1.add_argument('-o', '--output-file', default='-', help='A BED file containing fragment length summary statistics (mean, median, st. dev, min, max) over the intervals specified in the interval file.')
    parser_command3_1.add_argument('-q', '--quality-threshold', default=30, type=int, help='Minimum mapping quality threshold.')
    parser_command3_1.add_argument('-w', '--workers', default=1, type=int, help='Number of worker processes.')
    parser_command3_1.add_argument('-v', '--verbose', default=0, action='count', help='Enable verbose mode to display detailed processing information.')
    parser_command3_1.set_defaults(func=frag_length_intervals)

    parser_command13 = subparsers.add_parser('cleavage-profile', prog='finaletoolkit-cleavage-profile', description='Calculates cleavage proportion over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.')
    parser_command13.add_argument('input_file', help='Path to a BAM/SAM/CRAM/Fragment file containing fragment data.')
    parser_command13.add_argument('interval_file', help='Path to a BED file containing intervals to calculates cleavage proportion over.')
    parser_command13.add_argument('-o', '--output_file', default='-', help='A bigWig file containing the cleavage proportion results over the intervals specified in interval file.',)
    parser_command13.add_argument('-lo', '--fraction_low', default=120, type=int, help="Minimum length for a fragment to be included in cleavage proportion calculation.")
    parser_command13.add_argument('-hi', '--fraction_high', default=180, type=int, help="Maximum length for a fragment to be included in cleavage proportion calculation.")
    parser_command13.add_argument('-q', '--quality-threshold', default=20, type=int, help='Minimum mapping quality threshold.')
    parser_command13.add_argument('-l', '--left', default=0, type=int, help='Number of base pairs to subtract from start coordinate to create interval. Useful when dealing with BED files with only CpG coordinates.')
    parser_command13.add_argument('-r','--right', default=0, type=int, help='Number of base pairs to add to stop coordinate to create interval. Useful when dealing with BED files with only CpG coordinates.')
    parser_command13.add_argument('-w', '--workers', default=1, type=int, help='Number of worker processes.')
    parser_command13.add_argument('-v', '--verbose', default=0, action='count', help='Enable verbose mode to display detailed processing information.')
    parser_command13.set_defaults(func=_cli_cleavage_profile)

    parser_command4 = subparsers.add_parser('wps', prog='finaletoolkit-wps', description='Calculates Windowed Protection Score (WPS) over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.')
    parser_command4.add_argument('input_file', help='Path to a BAM/SAM/CRAM/Fragment file containing fragment data.')
    parser_command4.add_argument('site_bed', help='Path to a BED file containing intervals to calculate WPS over.')
    parser_command4.add_argument('-o', '--output_file', default='-', help='A bigWig file containing the WPS results over the intervals specified in interval file.')
    parser_command4.add_argument('-i', '--interval_size', default=5000, type=int, help='Size in bp of each interval in the interval file.')
    parser_command4.add_argument('-W', '--window_size', default=120, type=int, help='Size of the sliding window used to calculate WPS scores.')
    parser_command4.add_argument('-lo', '--fraction_low', default=120, type=int, help='Minimum length for a fragment to be included in WPS calculation.')
    parser_command4.add_argument('-hi', '--fraction_high', default=180, type=int,  help='Maximum length for a fragment to be included in WPS calculation.')
    parser_command4.add_argument('-q', '--quality_threshold', default=30, type=int, help="Minimum mapping quality threshold.")
    parser_command4.add_argument('-w', '--workers', default=1, type=int, help='Number of worker processes.')
    parser_command4.add_argument('-v', '--verbose', action='count', default=0, help='Enable verbose mode to display detailed processing information.')
    parser_command4.set_defaults(func=multi_wps)

    parser_command7 = subparsers.add_parser('adjust-wps', prog='finaletoolkit-adjust-wps', description='Adjusts raw Windowed Protection Score (WPS) by applying a median filter and Savitsky-Golay filter.')
    parser_command7.add_argument('input_file', help='A bigWig file containing the WPS results over the intervals specified in interval file.')
    parser_command7.add_argument('interval_file', help='Path to a BED file containing intervals to WPS was calculated over.')
    parser_command7.add_argument('genome_file', help='A .chrom.sizes file containing chromosome sizes.')
    parser_command7.add_argument('-o', '--output-file', default='-', help='A bigWig file containing the adjusted WPS results over the intervals specified in interval file.')
    parser_command7.add_argument('-m', '--median-window-size', default=1000, type=int, help='Size of the median filter window used to adjust WPS scores.')
    parser_command7.add_argument('-s', '--savgol-window-size', default=21, type=int, help='Size of the Savitsky-Golay filter window used to adjust WPS scores.')
    parser_command7.add_argument('-p', '--savgol-poly-deg', default=2, type=int, help='Degree polynomial for Savitsky-Golay filter.')
    parser_command7.add_argument('-w', '--workers', default=1, type=int, help='Number of worker processes.')
    parser_command7.add_argument('--mean', action='store_true', help='A mean filter is used instead of median.')
    parser_command7.add_argument('--subtract-edges', action='store_true', help='Take the median of the first and last 500 bases in a window and subtract from the whole interval.')
    parser_command7.add_argument('-v', '--verbose', action='count', help='Enable verbose mode to display detailed processing information.')
    parser_command7.set_defaults(func=adjust_wps)

    parser_command5 = subparsers.add_parser('delfi', prog='finaletoolkit-delfi', description='Calculates DELFI featrues over genome, returning information about (GC-corrected) short fragments, long fragments, DELFI ratio, and total fragments. NOTE: Due to some ad hoc implementation details, currently the only accepted reference genome is hg19.')
    parser_command5.add_argument('input_file', help="Path to a BAM/SAM/CRAM/Fragment file containing fragment data.")
    parser_command5.add_argument('autosomes', help="Tab-delimited file containing (1) autosome name and (2) integer length of chromosome in base pairs.")
    parser_command5.add_argument('reference_file', help="The .2bit file for the associate reference genome sequence used during alignment.")
    parser_command5.add_argument('bins_file', help="A BED file containing bins over which to calculate DELFI. To replicate Cristiano et al.'s methodology, use 100kb bins over human autosomes.")
    parser_command5.add_argument('-b', '--blacklist_file', help="BED file containing regions to ignore when calculating DELFI.")
    parser_command5.add_argument('-g', '--gap_file', help='BED4 format file containing columns for "chrom", "start", "stop", and "type". The "type" column should denote whether the entry corresponds to a "centromere", "telomere", or "short arm", and entries not falling into these categories are ignored. This information corresponds to the "gap" track for hg19 in the UCSC Genome Browser.')
    parser_command5.add_argument('-o', '--output_file', default='-', help='BED, bed.gz, TSV, or CSV file to write DELFI data to. If "-", writes to stdout.')
    parser_command5.add_argument('-W', '--window_size', default=5000000, type=int, help="Currently unused.")
    parser_command5.add_argument('-gc', '--gc_correct', action='store_true', help="Indicates whether to apply GC correction.")
    parser_command5.add_argument('-m', '--merge_bins', action='store_true', help="Indicates whether to merge bins to 5Mb size.")
    parser_command5.add_argument('-q', '--quality_threshold', default=30, type=int, help="Minimum mapping quality threshold.")
    parser_command5.add_argument('-w', '--workers', default=1, type=int, help="Number of worker processes.")
    parser_command5.add_argument('-v', '--verbose', action='count', default=0, help='Enable verbose mode to display detailed processing information.')
    parser_command5.set_defaults(func=delfi)

    parser_command9 = subparsers.add_parser('delfi-gc-correct', prog='finaletoolkit-delfi-gc-correct', description='Performs gc-correction on raw delfi data.')
    parser_command9.add_argument('input_file', help='BED file containing raw DELFI data. Raw DELFI data should only have columns for "contig", "start", "stop", "arm", "short", "long", "gc", "num_frags", "ratio".')
    parser_command9.add_argument('-o', '--output-file', default='-', help='BED to print GC-corrected DELFI fractions. If "-", will write to stdout.')
    parser_command9.add_argument('--header-lines', default=1, type=int, help='Number of header lines in BED.')
    parser_command9.add_argument('-v', '--verbose', action='count', help='Enable verbose mode to display detailed processing information.')
    parser_command9.set_defaults(func=cli_delfi_gc_correct)

    parser_command10 = subparsers.add_parser('end-motifs', prog='finaletoolkit-end-motifs', description="Measures frequency of k-mer 5' end motifs.")
    parser_command10.add_argument('input_file', help='Path to a BAM/SAM/CRAM/Fragment file containing fragment data.')
    parser_command10.add_argument('refseq_file', help='The .2bit file for the associate reference genome sequence used during alignment.')
    parser_command10.add_argument('-k', default=4, type=int, help='Length of k-mer.')
    parser_command10.add_argument('-o', '--output-file', default='-', help='TSV to print k-mer frequencies. If "-", will write to stdout.')
    parser_command10.add_argument('-q', '--quality-threshold', default=20, type=int, help='Minimum mapping quality threshold.')
    parser_command10.add_argument('-w', '--workers', default=1, type=int, help='Number of worker processes.')
    parser_command10.add_argument('-v', '--verbose', default=0, action='count', help='Enable verbose mode to display detailed processing information.')
    parser_command10.set_defaults(func=end_motifs)


    parser_command10a = subparsers.add_parser('interval-end-motifs', prog='finaletoolkit-interval-end-motifs', description="Measures frequency of k-mer 5' end motifs in each region specified in a BED file and writes data into a table.")
    parser_command10a.add_argument('input_file', help='Path to a BAM/SAM/CRAM/Fragment file containing fragment data.')
    parser_command10a.add_argument('refseq_file', help='The .2bit file for the associate reference genome sequence used during alignment.')
    parser_command10a.add_argument('intervals', help='Path to a BED file containing intervals to retrieve end motif frequencies over.')
    parser_command10a.add_argument('-k', default=4, type=int, help='Length of k-mer.')
    parser_command10a.add_argument('-lo', '--fraction-low', default=10, type=int, help='Minimum length for a fragment to be included in end motif frequency.')
    parser_command10a.add_argument('-hi', '--fraction-high', default=600, type=int, help='Maximum length for a fragment to be included in end motif frequency.')
    parser_command10a.add_argument('-o', '--output-file', default='-', help="Path to TSV or CSV file to write end motif frequencies to.")
    parser_command10a.add_argument('-q', '--quality-threshold', default=20, type=int, help='Minimum mapping quality threshold.')
    parser_command10a.add_argument('-w', '--workers', default=1, type=int, help='Number of worker processes.')
    parser_command10a.add_argument('-v', '--verbose', default=0, action='count', help='Enable verbose mode to display detailed processing information.')
    parser_command10a.set_defaults(func=interval_end_motifs)

    parser_command11 = subparsers.add_parser('mds', prog='finaletoolkit-mds', description='Reads k-mer frequencies from a file and calculates a motif diversity score (MDS) using normalized Shannon entropy as described by Jiang et al (2020).')
    parser_command11.add_argument('file_path', nargs='?', default='-', help='Tab-delimited or similar file containing one column for all k-mers a one column for frequency. Reads from stdin by default.')
    parser_command11.add_argument('-s', '--sep', default='\t', help='Separator used in tabular file.')
    parser_command11.add_argument('--header', default=0, type=int, help='Number of header rows to ignore. Default is 0')
    parser_command11.set_defaults(func=_cli_mds)

    parser_command11a = subparsers.add_parser('interval-mds', prog='finaletoolkit-interval-mds', description='Reads k-mer frequencies from a file and calculates a motif diversity score (MDS) for each interval using normalized Shannon entropy as described by Jiang et al (2020).')
    parser_command11a.add_argument('file_path', nargs='?', default='-', help='Tab-delimited or similar file containing one column for all k-mers a one column for frequency. Reads from stdin by default.')
    parser_command11a.add_argument('-s', '--sep', default='\t', help='Separator used in tabular file.')
    parser_command11a.add_argument('file_out', default='-', help='Path to the output BED/BEDGraph file containing MDS for each interval.')
    parser_command11a.set_defaults(func=_cli_interval_mds)

    parser_command6 = subparsers.add_parser('filter-bam', prog='finaletoolkit-filter-bam', description='Filters a BAM file so that all reads are in mapped pairs, exceed a certain MAPQ, are not flagged for quality, are read1, are not secondary or supplementary alignments, and are on the same reference sequence as the mate.')
    parser_command6.add_argument('input_file', help='Path to BAM file.')
    parser_command6.add_argument('-r', '--region-file', default=None, help='Only output alignments overlapping the intervals in this BED file will be included.')
    parser_command6.add_argument('-o', '--output-file', default='-', help='Output BAM file path.')
    parser_command6.add_argument('-q', '--quality_threshold', type=int, default=30, help='Minimum mapping quality threshold.')
    parser_command6.add_argument('-hi', '--fraction-high', type=int, default=None, help='Maximum length for a fragment to be included in output BAM.')
    parser_command6.add_argument('-lo', '--fraction-low', type=int, default=None, help='Minimum length for a fragment to be included in output BAM.')
    parser_command6.add_argument('-w', '--workers', type=int, default=1, help='Number of worker processes.')
    parser_command6.add_argument('-v', '--verbose', action='count', help='Enable verbose mode to display detailed processing information.')
    parser_command6.set_defaults(func=filter_bam)

    parser_command8 = subparsers.add_parser('agg-bw', prog='finaletoolkit-agg-wps', description='Aggregates a bigWig signal over constant-length intervals defined in a BED file.')
    parser_command8.add_argument('input_file', help=' A bigWig file containing signals over the intervals specified in interval file.')
    parser_command8.add_argument('interval_file', help='Path to a BED file containing intervals over which signals were calculated over.')
    parser_command8.add_argument('-o', '--output-file', default='-', help='A wiggle file containing the aggregate signal over the intervals specified in interval file.')
    parser_command8.add_argument('-m', '--median-window-size', default=0, type=int, help='Size of the median filter window used to adjust WPS scores. Only modify if aggregating WPS signals.')
    parser_command8.add_argument('-v', '--verbose', action='count', help='Enable verbose mode to display detailed processing information.')
    parser_command8.set_defaults(func=agg_bw)

    parser_command12 = subparsers.add_parser('gap-bed', prog='finaletoolkit-gap-bed', description='Creates a BED4 file containing centromeres, telomeres, and short-arm intervals, similar to the gaps annotation track for hg19 found on the UCSC Genome Browser (Kent et al 2002). Currently only supports hg19, b37, human_g1k_v37, hg38, and GRCh38', epilog='Gap is used liberally in this command, and in the case hg38/GRCh38, may refer to regions where there no longer are gaps in the reference sequence.')
    parser_command12.add_argument('reference_genome', choices=['hg19', 'b37','human_g1k_v37', 'hg38', 'GRCh38'], help='Reference genome to provide gaps for.')
    parser_command12.add_argument('output_file', help='Path to write BED file to. If "-" used, writes to stdout.')
    parser_command12.set_defaults(func=_cli_gap_bed)

 
    return parser



def main_cli():
    parser = main_cli_parser()

    args = parser.parse_args()
    try:
        function = args.func
        funcargs = vars(args)
        print(funcargs)
        funcargs.pop('func')

        function(**funcargs)
    except AttributeError:
        parser.print_help()
    
if __name__ == '__main__':
    main_cli()