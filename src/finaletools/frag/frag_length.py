from __future__ import annotations
import time
from typing import Union
from sys import stdout, stderr

import numpy as np
import pysam
from tqdm import tqdm

from finaletools.utils import not_read1_or_low_quality


def single_frag_length(input_file: Union[str, pysam.AlignedSegment],
                contig: str=None,
                start: int=None,
                stop: int=None,
                output_file: str=None, workers: int=1,
                quality_threshold: int=30,
                verbose: bool=False
                ) -> np.ndarray:
    """
    Return `np.ndarray` containing lengths of fragments in `input_file`
    that are above the quality threshold and are proper-paired reads.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM, SAM, or CRAM file containing paired-end fragment reads or
        its path. `AlignmentFile` must be opened in read mode.
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
        stderr.write("Finding frag lengths.\n")

    lengths = []    # list of fragment lengths
    if (type(input_file) == pysam.AlignmentFile):
        sam_file = input_file
        if (verbose):
            stderr.write('Counting reads\n')
        count = sam_file.count(contig=contig,
                               start=start,
                               stop=stop) if verbose else None
        if (verbose):
            stderr.write(f'{count} reads counted\n')
        # Iterating on each read in file in specified contig/chromosome
        for read1 in (tqdm(sam_file.fetch(contig=contig,
                                          start=start,
                                          stop=stop), total=count)
                      if verbose
                      else sam_file.fetch(contig=contig,
                                          start=start,
                                          stop=stop)):
            # Only select forward strand and filter out non-paired-end
            # reads and low-quality reads
            if (not_read1_or_low_quality(read1, quality_threshold)):
                pass
            else:
                # append length of fragment to list
                lengths.append(abs(read1.template_length))
    else:
        if (verbose):
            stderr.write(f'Opening {input_file}\n')
        with pysam.AlignmentFile(input_file) as sam_file:   # Import
            if (verbose):
                stderr.write('Counting reads')
            count = sam_file.count(contig=contig,
                                start=start,
                                stop=stop) if verbose else None
            if (verbose):
                stderr.write(f'{count} reads counted\n')
            # Iterating on each read in file in specified contig/chromosome
            for read1 in (tqdm(sam_file.fetch(contig=contig,
                                            start=start,
                                            stop=stop), total=count)
                        if verbose
                        else sam_file.fetch(contig=contig,
                                            start=start,
                                            stop=stop)):
                # Only select forward strand and filter out non-paired-end
                # reads and low-quality reads
                if (not_read1_or_low_quality(read1, quality_threshold)):
                    pass
                else:
                    # append length of fragment to list
                    lengths.append(abs(read1.template_length))

    # convert to array
    lengths = np.array(lengths)

    # check if output specified
    if (type(output_file) == str):
        if output_file.endswith(".bin"): # binary file
            with open(output_file, 'wt') as out:
                lengths.tofile(out)
        elif output_file == '-':
            for line in lengths:
                stdout.write(f'{line}\n')

        else:   # unaccepted file type
            raise ValueError(
                'output_file can only have suffixes .wig or .wig.gz.'
                )

    elif (output_file is not None):
        raise TypeError(
            f'output_file is unsupported type "{type(input_file)}". '
            'output_file should be a string specifying the path of the file '
            'to write output scores to.'
            )

    if (verbose):
        end_time = time.time()
        stderr.write(
            f'frag_length took {end_time - start_time} s to complete\n'
        )

    return lengths

def frag_length(

):
    return None