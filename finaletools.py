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

def low_quality_read_pairs(read, min_mapq=15): # min_mapq is synonymous to quality_threshold, copied from https://github.com/epifluidlab/cofragr/blob/master/python/frag_summary_in_intervals.py
    """
    Return `True` if the sequenced read described in `read` is not a mapped properly paired read with a Phred score exceeding `min_mapq`.

    Parameters
    ----------
    read : pysam.AlignedSegment
        Sequenced read to check for quality, perfect pairing and if it is mapped.
    min_mapq : int, optional
        Minimum Phred score for map quality of read. Defaults to 15.

    Returns
    -------
    is_low_quality : bool
        True if read is low quality, unmapped, not properly paired.
    """

    return read.is_unmapped or read.is_secondary or (not read.is_paired) \
           or read.mate_is_unmapped or read.is_duplicate or read.mapping_quality < min_mapq \
           or read.is_qcfail or read.is_supplementary or (not read.is_proper_pair) \
           or read.reference_name != read.next_reference_name

def frag_length(input_file, contig=None, output_file=None, threads=1, quality_threshold=15, verbose=False):
    """
    Return `np.ndarray` containing lengths of fragments in `input_file` that are
    above the quality threshold and are proper-paired reads.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, SAM, or CRAM file containing paired-end fragment reads or its 
        path. `AlignmentFile` must be opened in read mode.
    contig : string, optional
    output_file : string, optional
    quality_threshold : int, optional
    verbose : bool, optional

    Returns
    -------
    lengths : numpy.ndarray
        `ndarray` of fragment lengths from file and contig if 
        specified.
    """
    if (verbose):
        start_time = time.time()

    lengths = []    # list of fragment lengths
    if (type(input_file) == pysam.AlignmentFile):
        sam_file = input_file
        for read1 in sam_file.fetch(contig=contig): # Iterating on each read in file in specified contig/chromosome
            if read1.is_read2 or low_quality_read_pairs(read1, quality_threshold):  # Only select forward strand and filter out non-paired-end reads and low-quality reads
                pass
            else:
                lengths.append(abs(read1.template_length))  # append length of fragment to list
    else:
        with pysam.AlignmentFile(input_file) as sam_file:   # Import
            for read1 in sam_file.fetch(contig=contig): # Iterating on each read in file in specified contig/chromosome
                if read1.is_read2 or low_quality_read_pairs(read1, quality_threshold):  # Only select forward strand and filter out non-paired-end reads and low-quality reads
                    pass
                else:
                    lengths.append(abs(read1.template_length))  # append length of fragment to list
    # TODO: when given output file, print to file

    if (verbose):
        end_time = time.time()
        print(f'frag_coverage took {end_time - start_time} s to complete')
    return np.array(lengths)

# TODO: Read about pile-up

def frag_center_coverage(input_file, contig, start, stop, output_file=None, quality_threshold=15, verbose=False):
    """
    Return estimated fragment coverage over specified `contig` and region of 
    `input_file`. Uses an algorithm where the midpoints of fragments are calculated
    and coverage is tabulated from the midpoints that fall into the specified
    region. Not suitable for fragments of size approaching region size.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, SAM, or CRAM file containing paired-end fragment reads or its 
        path. `AlignmentFile` must be opened in read mode.
    contig : string
    start : int
    stop : int
    output_file : string, optional
    quality_threshold : int, optional
    verbose : bool, optional

    Returns
    -------
    coverage : int
        Estimated fragment coverage over contig and region.
    """
    # TODO: determine if reference (like as found in pysam) is necessary
    # TODO: consider including region string (like in pysam)

    if (verbose):
        start_time = time.time()

    coverage = 0 # initializing variable for coverage tuple outside of with statement

    if (type(input_file) == pysam.AlignmentFile):
        sam_file = input_file
        for read1 in sam_file.fetch(contig=contig): # Iterating on each read in file in specified contig/chromosome
            if read1.is_read2 or low_quality_read_pairs(read1, quality_threshold):  # Only select forward strand and filter out non-paired-end reads and low-quality reads
                pass
            else:
                # calculate mid-point of fragment
                center = read1.reference_start + read1.template_length // 2
                if ((center >= start) and (center < stop)):
                    coverage += 1
    else:
        with pysam.AlignmentFile(input_file, 'r') as sam_file:   # Import
            for read1 in sam_file.fetch(contig=contig): # Iterating on each read in file in specified contig/chromosome
                if read1.is_read2 or low_quality_read_pairs(read1, quality_threshold):  # Only select forward strand and filter out non-paired-end reads and low-quality reads
                    pass
                else:
                    # calculate mid-point of fragment
                    center = read1.reference_start + read1.template_length // 2
                    if ((center >= start) and (center < stop)):
                        coverage += 1
    if (verbose):
        end_time = time.time()
        print(f'frag_coverage took {end_time - start_time} s to complete')

    return coverage

def single_wps(alignment_file, contig, location, window_size=120, quality_threshold=15, verbose=False):
    """
    Return Windowed Protection Scores as specified in Snyder et al (2016).

    Parameters
    ----------
    input_file : pysam.AlignmentFile
        BAM, SAM, or CRAM file containing paired-end fragment reads.
    contig : str
        Contig that bp is on.
    location : int
        Location of bp at center of window.
    window_size : int, optional
        Size of window to calculate WPS. Defaults to k=120, corresponding to
        L-WPS.
    quality_threshold : int, optional
        Minimum read mapping quality. Defaults to 15.
    verbose : bool, optional

    Returns
    -------
    score : int
    """
    pass

def wps(input_file, contig, start, stop, output_file=None, window_size=120, quality_threshold=15, verbose=False):
    """
    Return Windowed Protection Scores as specified in Snyder et al (2016) over a
    region [start,stop).

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, SAM, or CRAM file containing paired-end fragment reads or its 
        path. `AlignmentFile` must be opened in read mode.
    contig : str
    start : int
    stop : int
    output_file : string, optional
    window_size : int, optional
    quality_threshold : int, optional
    verbose : bool, optional

    Returns
    -------
    scores : numpy.ndarray
    """
    # TODO: implement wps

    return None

# TODO: look through argparse args and fix them all
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
    parser_command3.add_argument('--region')   # samtools region string


    args = parser.parse_args()
    print(frag_center_coverage(**vars(args)))


