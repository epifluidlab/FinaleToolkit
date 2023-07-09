from __future__ import annotations
import tempfile as tf
import subprocess

import pysam



def filter_bam(
        input_file: str,
        region_file: str=None,
        output_file: str=None,
        quality_threshold: int=30,
        workers: int=1,
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
    workers : int, optional
    verbose : bool, optional

    Returns
    -------
    output_file : str
    """

    # create tempfile to contain filtered bam
    if output_file is None:
        _, output_file = tf.mkstemp(suffix='.bam')
    elif output_file.endswith('bam'):
        pass
    else:
        raise ValueError('output_file should have suffix .bam')

    try:
        # create temp dir to store intermediate sorted file
        temp_dir = tf.TemporaryDirectory()

        flag_filtered_bam = temp_dir.name + '/flag_filtered.bam'

        process = subprocess.run(
            f'samtools view {input_file} -F 3852 -f 66 -b -h -o '
            f'{flag_filtered_bam} -q {quality_threshold} -@ {workers}',
            shell=True,
            check=True
        ) 

        # filter for reads on different reference
        with pysam.AlignmentFile(flag_filtered_bam, 'rb') as in_file:
            with pysam.AlignmentFile(
                output_file,
                'wb',
                template=in_file
            ) as out_file:
                for read in in_file:
                    if read.reference_name == read.next_reference_name:
                        out_file.write(read)

    finally:
        temp_dir.cleanup()