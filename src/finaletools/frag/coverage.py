from __future__ import annotations
import time

import pysam
from numba import jit
from finaletools.utils import not_read1_or_low_quality


def coverage(input_file,
                         contig,
                         start=0,
                         stop=None,
                         output_file=None,
                         quality_threshold=30,
                         verbose=False):
    """
    Return estimated fragment coverage over specified `contig` and
    region of`input_file`. Uses an algorithm where the midpoints of
    fragments are calculated and coverage is tabulated from the
    midpoints that fall into the specified region. Not suitable for
    fragments of size approaching region size.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, SAM, or CRAM file containing paired-end fragment reads or
        its path. `AlignmentFile` must be opened in read mode.
    contig : string
    start : int
    stop : int
    output_file : string, optional
    quality_threshold : int, optional
    verbose : bool, optional

    Returns
    -------
    coverage : int
        Fragment coverage over contig and region.
    """
    # TODO: determine if reference (like as found in pysam) is necessary
    # TODO: consider including region string (like in pysam)

    if (verbose):
        start_time = time.time()

    # initializing variable for coverage tuple outside of with statement
    coverage = 0

    if (type(input_file) == pysam.AlignmentFile):
        sam_file = input_file
        # Iterating on each read in file in specified contig/chromosome
        for read1 in sam_file.fetch(
            contig=contig,
            start=None if (tempstart:=start-500) < 0 else tempstart,
            stop=stop+500 if stop is not None else None):
            # Only select forward strand and filter out non-paired-end
            # reads and low-quality reads
            if (not_read1_or_low_quality(read1, quality_threshold)):
                pass
            else:
                # calculate mid-point of fragment
                center = read1.reference_start + read1.template_length // 2
                if start is None and stop is None:
                    coverage += 1
                elif start is None:
                    if center < stop:
                        coverage += 1
                elif stop is None:
                    if center >= start:
                        coverage += 1
                elif (center >= start) and (center < stop):
                    coverage += 1
                else:
                    pass
    else:
        with pysam.AlignmentFile(input_file, 'r') as sam_file:
            # Iterating on each read in file in
            # specified contig/chromosome
            for read1 in sam_file.fetch(
                contig=contig,
                start=0 if (tempstart:=start-500) < 0 else tempstart,
                stop=stop+500 if stop is not None else None):
                # Only select forward strand and filter out
                # non-paired-end reads and low-quality reads
                if not_read1_or_low_quality(read1, quality_threshold):
                    pass
                else:
                    # calculate mid-point of fragment
                    center = read1.reference_start + read1.template_length // 2
                    if start is None and stop is None:
                        coverage += 1
                    elif start is None:
                        if center < stop:
                            coverage += 1
                    elif stop is None:
                        if center >= start:
                            coverage += 1
                    elif (center >= start) and (center < stop):
                        coverage += 1
                    else:
                        pass
    if (verbose):
        end_time = time.time()
        print(f'frag_coverage took {end_time - start_time} s to complete')

    return coverage