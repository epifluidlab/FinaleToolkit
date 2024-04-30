from __future__ import annotations
import tempfile as tf
import subprocess
import traceback
from sys import stderr
from os.path import isdir
from shutil import rmtree

import pysam



def filter_bam(
        input_file: str,
        region_file: str=None,
        output_file: str=None,
        fraction_high: int=None,
        fraction_low: int=None,
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
    elif output_file.endswith('bam') or output_file == '-':
        pass
    else:
        raise ValueError('output_file should have suffix .bam')


        # create temp dir to store intermediate sorted file
    try:
        temp_dir = tf.TemporaryDirectory()

        flag_filtered_bam = temp_dir.name + '/flag_filtered.bam'

        samtools_command = (
            f'samtools view {input_file} -F 3852 -f 66 -b -h -o '
            f'{flag_filtered_bam} -q {quality_threshold} -@ {workers}'
        )

        if region_file is not None:
            samtools_command += f' -M -L {region_file}'

        try:
            process1 = subprocess.run(samtools_command, shell=True, check=True)
        except Exception as e:
            traceback.print_exc()
            exit(1)

        # supress index file warning
        save = pysam.set_verbosity(0)

        # filter for reads on different reference
        with pysam.AlignmentFile(flag_filtered_bam, 'rb') as in_file:
            with pysam.AlignmentFile(
                output_file,
                'wb',
                template=in_file
            ) as out_file:
                if fraction_high is None and fraction_low is None:
                    for read in in_file:
                        if read.reference_name == read.next_reference_name:
                            out_file.write(read)
                elif fraction_high is None:
                    for read in in_file:
                        if (read.reference_name == read.next_reference_name
                            and read.template_length >= fraction_low
                        ):
                            out_file.write(read)
                elif fraction_low is None:
                    for read in in_file:
                        if (read.reference_name == read.next_reference_name
                            and read.template_length <= fraction_high
                        ):
                            out_file.write(read)
                else:
                    for read in in_file:
                        if (read.reference_name == read.next_reference_name
                            and read.template_length >= fraction_low
                            and read.template_length <= fraction_high
                        ):
                            out_file.write(read)
    finally:
        temp_dir.cleanup()
        pysam.set_verbosity(save)

    if output_file != '-':
        # generate index for output_file
        try:
            process3 = subprocess.run(
                f'samtools index {output_file} {output_file}.bai',
                shell=True,
                check=True
            )
        except Exception as e:
            traceback.print_exc()
            exit(1)