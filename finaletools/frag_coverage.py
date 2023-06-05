"""
Author: James Li
Created: 6/2/23
PI: Yaping Liu

Description:
Python script to calculate coverage of a window given a BAM file.

"""

import pysam
import argparse
import gzip
import numpy as np
from multiprocessing.pool import Pool
import time

# TODO: Read about pile-up
# TODO: finish frag coverage
def frag_coverage(input_file, contig, output_file=None, reference=None, start=None, stop=None, end=None, region=None, quality_threshold=15, read_callback='all', verbose=False):
    coverage = None # initializing variable for coverage tuple outside of with statement
    with pysam.AlignmentFile(input_file, 'r') as sam_file:
        False
    return coverage

if __name__ == '__main__':
    p = argparse.ArgumentParser(prog='frag_coverage',
                    description='Calculates coverage of a region given a BAM file',
                    epilog='Bottom text') # TODO: fix this bottom text
    p.add_argument('contig')   # chromosome of window
    p.add_argument('input_file')    # input bam file to calculate coverage from
    p.add_argument('--output_file') # optional output text file to print coverage in
    p.add_argument('--reference')   # synonymous to contig
    p.add_argument('--start', type=int)   # inclusive location of region start in 0-based coordinate system. If not included, will start at the beginning of the chromosome
    p.add_argument('--stop', type=int)   # exclusive location of region end in 0-based coordinate system. If not included, will end at the end of the chromosome
    p.add_argument('--end', type=int)   # synonymous to stop
    p.add_argument('--region')   # samtools region string
    p.add_argument('--quality_threshold', default=15, type=int)   # minimum phred score for a base to be counted TODO: implement threshold
    p.add_argument('--read_callback', default='all')    # TODO: implement read callback
    p.add_argument('-v', '--verbose', default=False, type=bool)    # TODO: add verbose mode
    args = p.parse_args()
    print(frag_coverage(**vars(args)))


