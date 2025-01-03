from __future__ import annotations
import tempfile as tf
import subprocess
import traceback
import warnings
import gzip
import pysam

def filter_file(
        input_file: str,
        whitelist_file: str | None = None,
        blacklist_file: str | None = None,
        output_file: str | None = None,
        min_length: int | None = None,
        max_length: int | None = None,
        quality_threshold: int = 30,
        workers: int = 1,
        verbose: bool = False,
        fraction_low: int | None = None,
        fraction_high: int | None = None,):
    """
    Accepts the path to a BAM, CRAM, or BED file and creates a filtered version.

    Filter reads/intervals based on exceeding the specified quality threshold,
    intersections with a region in the region bed (if provided), and read length.

    For BAM/CRAM files, it also filters reads based on being read1 in a proper pair.

    Parameters
    ----------
    input_file : str
        Path string to the input BAM, CRAM, or BED file.
    whitelist_file : str, optional
        Path to a BED file defining regions to include.
    blacklist_file : str, optional
        Path to a BED file defining regions to exclude.
    output_file : str, optional
        Path to the output filtered file. If None, a temporary file is created.
    min_length : int, optional
        Minimum length for reads/intervals
    max_length : int, optional
        Maximum length for reads/intervals
    quality_threshold : int, optional
        Minimum mapping quality score
    workers : int, optional
        Number of worker threads for samtools.
    verbose : bool, optional
        Default is False
    fraction_low : int, optional
        Deprecated alias for min_length
    fraction_high : int, optional
        Deprecated alias for max_length

    Returns
    -------
    output_file : str
        Path to the filtered output file.
    """

    if verbose:
        print(
        f"""
        input_file: {input_file}
        whitelist_file: {whitelist_file}
        blacklist_file: {blacklist_file}
        output_file: {output_file}
        min_length: {min_length}
        max_length: {max_length}
        quality_threshold: {quality_threshold}
        workers: {workers}
        verbose: {verbose}     
        """
        )
        
    # Pass aliases and check for conflicts
    if fraction_low is not None and min_length is None:
        min_length = fraction_low
        warnings.warn("fraction_low is deprecated. Use min_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
    elif fraction_low is not None and min_length is not None:
        warnings.warn("fraction_low is deprecated. Use min_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
        raise ValueError(
            'fraction_low and min_length cannot both be specified')

    if fraction_high is not None and max_length is None:
        max_length = fraction_high
        warnings.warn("fraction_high is deprecated. Use max_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
    elif fraction_high is not None and max_length is not None:
        warnings.warn("fraction_high is deprecated. Use max_length instead.",
                      category=DeprecationWarning,
                      stacklevel=2)
        raise ValueError(
            'fraction_high and max_length cannot both be specified.')

    if input_file.endswith(".gz"):
        suffix = ".gz"
    elif input_file.endswith(".bam"):
        suffix = ".bam"
    elif input_file.endswith(".cram"):
        suffix = ".cram"
    else:
        raise ValueError('Input file should have suffix .bam, .cram, or .gz')      

    # create tempfile to contain filtered output
    if output_file is None:
        _, output_file = tf.mkstemp(suffix=suffix)
    elif not output_file.endswith(suffix) and output_file != '-':
        raise ValueError('Output file should share same suffix as input file.')

    if input_file.endswith(('.bam', '.cram')):
        # create temp dir to store intermediate sorted file
        try:
            with tf.TemporaryDirectory() as temp_dir:
                flag_filtered_bam = f'{temp_dir}/flag_filtered{suffix}'
                samtools_command = (
                    f'samtools view {input_file} -F 3852 -f 3 -b -h -o '
                    f'{flag_filtered_bam} -q {quality_threshold} -@ {workers}'
                )

                if whitelist_file is not None:
                    samtools_command += f' -M -L {whitelist_file}'
                if blacklist_file is not None:
                    samtools_command += f' -U {blacklist_file}'
                try:
                    subprocess.run(samtools_command, shell=True, check=True)
                except Exception:
                    traceback.print_exc()
                    exit(1)

                # suppress index file warning
                save = pysam.set_verbosity(0)
                
                # filter for reads on different reference and length
                with pysam.AlignmentFile(flag_filtered_bam, 'rb',threads=workers//3) as in_file:
                    with pysam.AlignmentFile(
                        output_file, 'wb', template=in_file, threads=workers-workers//3) as out_file:
                        for read in in_file:
                            if (
                                read.reference_name == read.next_reference_name
                                and (max_length is None
                                     or read.template_length <= max_length)
                                and (min_length is None
                                     or read.template_length >= min_length)
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
    elif input_file.endswith('.gz'):
        with tf.TemporaryDirectory() as temp_dir:
            with gzip.open(input_file, 'r') as infile, open(f"{temp_dir}/output1.bed", 'w') as outfile:
                mapq_column = 0 # 1-index for sanity when comparing with len()
                for line in infile:
                    line = line.decode('utf-8')
                    parts = line.strip().split('\t')
                    if len(parts) < max(mapq_column,4) or line.startswith('#'):
                        continue

                    if mapq_column == 0:
                        if parts[4-1].isnumeric():
                            mapq_column = 4
                        elif len(parts) >= 5 and parts[5-1].isnumeric():
                            mapq_column = 5
                        else:
                            continue
                    try:
                        start = int(parts[1])
                        end = int(parts[2])
                        length = end - start
                        score = None
                        try:
                            score = float(parts[mapq_column-1])
                        except ValueError:
                            pass

                        passes_length_restriction = True
                        
                        if min_length is not None and length < min_length:
                            passes_length_restriction = False

                        if max_length is not None and length > max_length:
                            passes_length_restriction = False

                        passes_quality_restriction = True
                        if score is None or score < quality_threshold:
                            passes_quality_restriction = False

                        if passes_length_restriction and passes_quality_restriction:
                            outfile.write(line)
                    except ValueError:
                        continue
                if whitelist_file is not None:
                    try:
                        subprocess.run(
                            f"bedtools intersect -f 0.50 -a {temp_dir}/output1.bed -b {whitelist_file} > {temp_dir}/output2.bed", 
                            shell=True, 
                            check=True
                            )
                    except Exception:
                        traceback.print_exc()
                        exit(1)
                else:
                    subprocess.run(f"mv {temp_dir}/output1.bed {temp_dir}/output2.bed", shell=True, check=True)

                if blacklist_file is not None:
                    try:
                        subprocess.run(
                            f"bedtools intersect -v -f 0.50 -a {temp_dir}/output2.bed -b {blacklist_file} > {temp_dir}/output3.bed", 
                            shell=True, 
                            check=True
                            )
                    except Exception:
                        traceback.print_exc()
                        exit(1)  
                else:
                    subprocess.run(f"mv {temp_dir}/output2.bed {temp_dir}/output3.bed", shell=True, check=True)
                try:
                    subprocess.run(
                        f"bgzip -@ {workers} -c {temp_dir}/output3.bed > {output_file}", 
                        shell=True, 
                        check=True
                        )
                except Exception:
                    traceback.print_exc()
                    exit(1)    
                if output_file != '-':
                    # generate index for output_file
                    try:
                        subprocess.run(
                            f'tabix -p bed {output_file}',
                            shell=True,
                            check=True
                        )
                    except Exception:
                        traceback.print_exc()
                        exit(1)                        
    else:
        raise ValueError("Input file must be a BAM, CRAM, or BED file.")
    return output_file
