"""
Author: James Li
Created: 6/2/23
PI: Yaping Liu

Description:
Python script to calculate fragment features given a BAM file.

"""

import pysam
import argparse
import gzip
import numpy as np
from multiprocessing.pool import Pool
import time
from tempfile import TemporaryDirectory

def low_quality_read_pairs(read, min_mapq): # min_mapq is synonymous to quality_threshold
    return read.is_unmapped or read.is_secondary or (not read.is_paired) \
           or read.mate_is_unmapped or read.is_duplicate or read.mapping_quality < min_mapq \
           or read.is_qcfail or read.is_supplementary or (not read.is_proper_pair) \
           or read.reference_name != read.next_reference_name

def frag_length(input_file, contig=None, output_file=None, threads=1, quality_threshold=15, verbose=False):
    lengths = []    # list of fragment lengths
    with pysam.AlignmentFile(input_file) as sam_file:   # Import
        for read1 in sam_file.fetch(contig=contig): # Iterating on each read in file in specified contig/chromosome
            if read1.is_read2 or low_quality_read_pairs(read1, quality_threshold):  # Filter out non-paired-end reads and low-quality reads
                pass
            else:
                lengths.append(abs(read1.template_length))  # append length of fragment to list
    # TODO: when given output file, print to file
    return np.array(lengths)

# TODO: Read about pile-up

# TODO: finish frag coverage
def frag_coverage(input_file, contig, output_file=None, reference=None, start=None, stop=None, end=None, region=None, quality_threshold=15, read_callback='all', verbose=False):
    coverage = None # initializing variable for coverage tuple outside of with statement
    with pysam.AlignmentFile(input_file, 'r') as sam_file:
        False
    return coverage

def wps():
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='finaletools',
                    description='Calculates fragmentation features given a CRAM/BAM/SAM file',
                    epilog='')
    subparsers = parser.add_subparsers(title='subcommands', description="Subcommands of finaletools", dest='subcommand')

    # Common arguments
    parser.add_argument('contig')   # chromosome of window
    parser.add_argument('input_file')    # input bam file to calculate coverage from
    parser.add_argument('--output_file') # optional output text file to print coverage in
    parser.add_argument('--reference')   # synonymous to contig
    parser.add_argument('--quality_threshold', default=15, type=int)

    # Subcommand 1: frag-coverage
    parser_command1 = subparsers.add_parser(prog='frag-converage',
                                            description='Calculates fragmentation coverage over a region given a CRAM/BAM/SAM file')
    parser_command1.add_argument('--start', type=int)   # inclusive location of region start in 0-based coordinate system. If not included, will start at the beginning of the chromosome
    parser_command1.add_argument('--stop', type=int)   # exclusive location of region end in 0-based coordinate system. If not included, will end at the end of the chromosome
    parser_command1.add_argument('--end', type=int)   # synonymous to stop
    parser_command1.add_argument('--region')   # samtools region string
    parser_command1.add_argument('--method', default="frag-center")
    # parser_command1.add_argument('--quality_threshold', default=15, type=int)   # minimum phred score for a base to be counted TODO: implement threshold
    # parser_command1.add_argument('--read_callback', default='all')    # TODO: implement read callback
    parser_command1.add_argument('-v', '--verbose', default=False, type=bool)    # TODO: add verbose mode
    
    # Subcommand 2: frag-length
    parser_command2 = subparsers.add_parser(prog='frag-length',
                                            description='Calculates fragment lengths given a CRAM/BAM/SAM file')
    
    # Subcommand 2: wps()
    parser_command3 = subparsers.add_parser(prog='frag-length',
                                            description='Calculates Windowed Protection Score over a region given a CRAM/BAM/SAM file')
    parser_command3.add_argument('--start', type=int)   # inclusive location of region start in 0-based coordinate system. If not included, will start at the beginning of the chromosome
    parser_command3.add_argument('--stop', type=int)   # exclusive location of region end in 0-based coordinate system. If not included, will end at the end of the chromosome
    parser_command3.add_argument('--end', type=int)   # synonymous to stop
    parser_command3.add_argument('--region')   # samtools region string


    args = parser.parse_args()
    print(frag_coverage(**vars(args)))


