from __future__ import annotations
import subprocess
import traceback
from os.path import isdir
import pysam
import tempfile as tf

def filter_bam(
    input_file: str,
    region_file: str,
    output_file: str,
    max_length: int,
    min_length: int,
    quality_threshold: int,
    workers: int,
    verbose: bool
):

    if output_file is None:
        _, output_file = tf.mkstemp(suffix='.bam')
    elif not output_file.endswith('bam') and output_file != '-':
        raise ValueError('Output file should have suffix .bam')

    # create temp dir to store intermediate sorted file
    try:
        with tf.TemporaryDirectory() as temp_dir:
            flag_filtered_bam = f'{temp_dir}/flag_filtered.bam'
            samtools_command = (
                f'samtools view {input_file} -F 3852 -f 66 -b -h -o '
                f'{flag_filtered_bam} -q {quality_threshold} -@ {workers}'
            )

            if region_file is not None:
                samtools_command += f' -M -L {region_file}'

            try:
                subprocess.run(samtools_command, shell=True, check=True)
            except Exception:
                traceback.print_exc()
                exit(1)

            # suppress index file warning
            save = pysam.set_verbosity(0)

            # filter for reads on different reference
            with pysam.AlignmentFile(flag_filtered_bam, 'rb') as in_file:
                with pysam.AlignmentFile(output_file, 'wb', template=in_file) as out_file:
                    for read in in_file:
                        if (
                            read.reference_name == read.next_reference_name and
                            (max_length is None or read.template_length <= max_length) and
                            (min_length is None or read.template_length >= min_length)
                        ):
                            out_file.write(read)

    finally:
        pysam.set_verbosity(save)

    if output_file != '-':
        # generate index for output_file
        try:
            subprocess.run(
                f'samtools index {output_file} {output_file}.bai',
                shell=True,
                check=True
            )
        except Exception:
            traceback.print_exc()
            exit(1)
