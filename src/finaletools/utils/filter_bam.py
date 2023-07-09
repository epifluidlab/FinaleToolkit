from __future__ import annotations
import tempfile as tf

import pysam



def filter_bam(
        input_file: str,
        region_file: str=None,
        output_file: str=None,
        quality_threshold: int=30,
        verbose: bool=False):
    """
    Accepts the path to a BAM file and creates a bam file where all
    reads are read1 in a proper pair, exceed the specified quality
    threshold, do not intersect a region in the given blacklist
    file, and intersects with a region in the region bed.

    Parameters
    ----------
    input_bam : str
        Path string or AlignmentFile pointing to the BAM file to be
        filtered.
    region_file : str, option
    output_file : str, optional
    quality_threshold : int, optional
    verbose : bool, optional

    Returns
    -------
    output_file : str
    """

    # create tempfile to contain filtered BED
    if output_file is None:
        _, output_file = tf.mkstemp(suffix='.bed')
    elif output_file.endswith('bed'):
        pass
    else:
        raise ValueError('output_file should have suffix .bed')

    temp_dir = tf.mkdtemp()

    flag_filtered_bam = temp_dir + '/flag_filtered.bam'

    print(flag_filtered_bam)

    pysam.view('-b', '-h', '-o', flag_filtered_bam, '-q',
               str(quality_threshold), '-f', '2', input_file)

    with pysam.AlignmentFile(flag_filtered_bam, 'rb') as file:
        for line in file.head(5):
            pass
