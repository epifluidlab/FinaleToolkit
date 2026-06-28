"""
Filter a BAM/CRAM or tabix BED fragment file by quality, length, and region.

BAM/CRAM filtering shells out to bedtools/samtools; tabix BED filtering streams
lines and shells out to bedtools/bgzip/tabix.
"""
from __future__ import annotations

import gzip
import logging
import subprocess
import tempfile as tf
import warnings

import pysam

__all__ = ["filter_file"]

_VALID_SUFFIXES = (".gz", ".bam", ".cram", ".bgz")


def validate_deprecated_args(old_arg, new_arg, old_name, new_name):
    """Resolve a deprecated alias against its replacement argument."""
    if old_arg is not None:
        warnings.warn(
            f"{old_name} is deprecated. Use {new_name} instead.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        if new_arg is None:
            return old_arg
        raise ValueError(
            f"{old_name} and {new_name} cannot both be specified."
        )
    return new_arg


def validate_input_file(input_file: str) -> str:
    """Return the recognized suffix of ``input_file`` or raise ``ValueError``."""
    if not any(input_file.endswith(suffix) for suffix in _VALID_SUFFIXES):
        raise ValueError(
            "Input file should have one of the following suffixes: "
            f"{', '.join(_VALID_SUFFIXES)}"
        )
    return next(suffix for suffix in _VALID_SUFFIXES if input_file.endswith(suffix))


def run_subprocess(
    cmd: str,
    error_msg: str = "Command failed",
    verbose: bool = False,
    logger=None,
) -> None:
    """Run a shell command, logging and re-raising on failure."""
    try:
        if verbose:
            logger.info(f"Running: {cmd}")
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"{error_msg}: {str(e)}")
        raise


def filter_bed_entries(
    infile, min_length=None, max_length=None, quality_threshold=30
):
    """Yield BED lines passing the length and MAPQ filters."""

    def get_mapq_col(length):
        if length < 5:
            raise ValueError(
                "There are not enough columns in the BED file to determine "
                "the MAPQ column"
            )
        return 3 if length == 5 else 4

    for line in infile:
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")
        mapq_column = get_mapq_col(len(parts))

        try:
            start = int(parts[1])
            end = int(parts[2])
            length = end - start
            score = float(parts[mapq_column])

            if (
                (min_length is None or length >= min_length)
                and (max_length is None or length <= max_length)
                and score >= quality_threshold
            ):
                yield line
        except (ValueError, IndexError):
            continue


def filter_file(
    input_file: str,
    whitelist_file: str | None = None,
    blacklist_file: str | None = None,
    output_file: str | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    intersect_policy: str = "midpoint",
    quality_threshold: int = 30,
    workers: int = 1,
    verbose: bool = False,
    fraction_low: int | None = None,
    fraction_high: int | None = None,
) -> str:
    """Create a filtered copy of a BAM/CRAM or tabix BED fragment file.

    Reads/intervals are kept when they exceed ``quality_threshold``, fall
    within the length bounds, and (for whitelist/blacklist BEDs) intersect per
    ``intersect_policy``.  BAM/CRAM additionally keeps only read1 of proper
    pairs.

    Parameters
    ----------
    input_file : str
        Input BAM/CRAM/`.gz`/`.bgz` path.
    whitelist_file : str, optional
        BED of regions to keep.
    blacklist_file : str, optional
        BED of regions to exclude.
    output_file : str, optional
        Output path (``"-"`` for stdout).
    min_length, max_length : int, optional
        Fragment-length bounds.
    intersect_policy : {"midpoint", "any"}, optional
        Whitelist/blacklist intersection policy (default ``"midpoint"``).
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    workers : int, optional
        Threads for samtools/bgzip (default 1).
    verbose : bool, optional
        Print configuration and enable logging.
    fraction_low, fraction_high : int, optional
        Deprecated aliases for ``min_length``/``max_length``.

    Returns
    -------
    str
        The output path.
    """
    logger = logging.getLogger(__name__)
    if verbose:
        print(
            f"""
        input_file: {input_file}
        whitelist_file: {whitelist_file}
        blacklist_file: {blacklist_file}
        output_file: {output_file}
        min_length: {min_length}
        max_length: {max_length}
        intersect_policy: {intersect_policy}
        quality_threshold: {quality_threshold}
        workers: {workers}
        verbose: {verbose}
        """
        )
        logging.basicConfig(level=logging.INFO)

    min_length = validate_deprecated_args(
        fraction_low, min_length, "fraction_low", "min_length"
    )
    max_length = validate_deprecated_args(
        fraction_high, max_length, "fraction_high", "max_length"
    )

    suffix = validate_input_file(input_file)

    if intersect_policy == "midpoint":
        intersect_param = "-f 0.500"
    elif intersect_policy == "any":
        intersect_param = ""
    else:
        raise ValueError("intersect_policy must be 'midpoint' or 'any'")

    pysam.set_verbosity(pysam.set_verbosity(0))

    with tf.TemporaryDirectory() as temp_dir:
        temp_1 = f"{temp_dir}/output1{suffix}"
        temp_2 = f"{temp_dir}/output2{suffix}"
        temp_3 = f"{temp_dir}/output3{suffix}"
        if input_file.endswith((".bam", ".cram")):
            _filter_alignment(
                input_file,
                whitelist_file,
                blacklist_file,
                output_file,
                min_length,
                max_length,
                intersect_policy,
                intersect_param,
                quality_threshold,
                workers,
                verbose,
                logger,
                temp_1,
                temp_2,
                temp_3,
            )
        elif input_file.endswith((".gz", ".bgz")):
            _filter_bed(
                input_file,
                whitelist_file,
                blacklist_file,
                output_file,
                min_length,
                max_length,
                intersect_policy,
                quality_threshold,
                workers,
                verbose,
                logger,
                temp_1,
                temp_2,
                temp_3,
            )
    return output_file


def _filter_alignment(
    input_file,
    whitelist_file,
    blacklist_file,
    output_file,
    min_length,
    max_length,
    intersect_policy,
    intersect_param,
    quality_threshold,
    workers,
    verbose,
    logger,
    temp_1,
    temp_2,
    temp_3,
) -> None:
    """Filter a BAM/CRAM via bedtools/samtools and a length pass."""
    if whitelist_file:
        run_subprocess(
            f"bedtools intersect -abam {input_file} -b {whitelist_file} "
            f"{intersect_param} > {temp_1} && samtools index {temp_1}",
            error_msg="Whitelist filtering failed",
            verbose=verbose,
            logger=logger,
        )
    else:
        run_subprocess(f"cp {input_file} {temp_1}", verbose=verbose, logger=logger)

    if blacklist_file:
        intersect_param = "-f 0.500" if intersect_policy == "midpoint" else ""
        run_subprocess(
            f"bedtools intersect -abam {temp_1} -b {blacklist_file} -v "
            f"{intersect_param} > {temp_2} && samtools index {temp_2}",
            error_msg="Blacklist filtering failed",
            verbose=verbose,
            logger=logger,
        )
    else:
        run_subprocess(f"mv {temp_1} {temp_2}", verbose=verbose, logger=logger)

    run_subprocess(
        f"samtools view {temp_2} -F 3852 -f 3 -b -h -o {temp_3} "
        f"-q {quality_threshold} -@ {workers}",
        error_msg="Quality filtering failed",
        verbose=verbose,
        logger=logger,
    )

    # Length filter (template length, same reference for both mates).
    pysam.set_verbosity(0)
    with pysam.AlignmentFile(temp_3, "rb", threads=workers // 3) as in_file:
        with pysam.AlignmentFile(
            output_file, "wb", template=in_file, threads=workers - workers // 3
        ) as out_file:
            for read in in_file:
                if (
                    read.reference_name == read.next_reference_name
                    and (max_length is None or read.template_length <= max_length)
                    and (min_length is None or read.template_length >= min_length)
                ):
                    out_file.write(read)

    if output_file != "-":
        run_subprocess(
            f"samtools index {output_file}",
            error_msg="Index creation failed",
            verbose=verbose,
            logger=logger,
        )


def _filter_bed(
    input_file,
    whitelist_file,
    blacklist_file,
    output_file,
    min_length,
    max_length,
    intersect_policy,
    quality_threshold,
    workers,
    verbose,
    logger,
    temp_1,
    temp_2,
    temp_3,
) -> None:
    """Filter a tabix BED fragment file via streaming + bedtools/bgzip/tabix."""
    with gzip.open(input_file, "rt") as infile, open(temp_1, "w") as outfile:
        for line in filter_bed_entries(
            infile, min_length, max_length, quality_threshold
        ):
            outfile.write(line)
            outfile.flush()

        if whitelist_file:
            intersect_param = "-f 0.500" if intersect_policy == "midpoint" else ""
            run_subprocess(
                f"bedtools intersect -a {temp_1} -b {whitelist_file} "
                f"{intersect_param} > {temp_2}",
                error_msg="Whitelist filtering failed",
                verbose=verbose,
                logger=logger,
            )
        else:
            run_subprocess(f"mv {temp_1} {temp_2}", verbose=verbose, logger=logger)

        if blacklist_file:
            intersect_param = "-f 0.500" if intersect_policy == "midpoint" else ""
            run_subprocess(
                f"bedtools intersect -v -a {temp_2} -b {blacklist_file} "
                f"{intersect_param} > {temp_3}",
                error_msg="Blacklist filtering failed",
                verbose=verbose,
                logger=logger,
            )
        else:
            run_subprocess(f"mv {temp_2} {temp_3}", verbose=verbose, logger=logger)

        run_subprocess(
            f"bgzip -@ {workers} -c {temp_3} > {output_file}",
            error_msg="Compression failed",
            verbose=verbose,
            logger=logger,
        )

        if output_file != "-":
            run_subprocess(
                f"tabix -p bed {output_file}",
                error_msg="Index creation failed",
                verbose=verbose,
                logger=logger,
            )
