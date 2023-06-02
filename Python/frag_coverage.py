"""
Author: James Li
Created: 6/2/23
PI: Yaping Liu

Description:
Python script to calculate coverage of a window or nucleotide given a BAM file.
"""

import pysam
import argparse
import gzip
import numpy as np
from multiprocessing.pool import Pool
import time

def frag_coverage(bam_file, chrom, start, end):
    
    return None

if __name__ == '__main__':
    p = argparse.ArgumentParser(prog='frag_coverage',
                    description='Calculates coverage of a window or nucleotide given a BAM file',
                    epilog='Bottom text') # TODO: fix this bottom text
    p.add_argument('--bam_file')    # input bam file to calculate coverage from
    p.add_argument('--chrom')   # chromosome of window
    p.add_argument('--start')   # location of start nt or single nt in 1-based coordinate system TODO: verify that pysam uses 1-based
    p.add_argument('--end')   # location of end nt. Set to 'start' for single nt
    args = p.parse_args()
    print(args.bam_file, args.chrom, args.start, args.end)


