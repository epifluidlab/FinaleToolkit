"""
Aggregate WPS over many sites in a BED file (Snyder et al., 2016).
"""
from __future__ import annotations

import gzip
import time
import warnings
from multiprocessing.pool import Pool
from os import PathLike
from pathlib import Path
from sys import stderr, stdin
from typing import Union

import numpy as np
import pyBigWig as pbw
import pysam

from finaletoolkit.frag._wps import wps
from finaletoolkit.utils import chrom_sizes_to_list
from finaletoolkit.utils.typing import ChromSizes, FragFile, Intervals

__all__ = ["multi_wps"]


def _wps_star(args):
    """Unpack a tuple of arguments and call :func:`wps` (for ``Pool.imap``)."""
    return wps(*args)


def multi_wps(
    input_file: FragFile,
    site_bed: Intervals,
    chrom_sizes: ChromSizes | None = None,
    output_file: str | None = None,
    window_size: int = 120,
    interval_size: int = 5000,
    min_length: int = 120,
    max_length: int = 180,
    quality_threshold: int = 30,
    workers: int = 1,
    verbose: Union[bool, int] = 0,
    fraction_low: int | None = None,
    fraction_high: int | None = None,
    reference_file: str | Path | None = None,
) -> str | None:
    """Aggregate WPS over sites in a BED file.

    Each BED interval is replaced by a ``interval_size``-bp window centered on
    its midpoint; WPS is computed per window and written to a bigWig (or
    bed.gz/bedGraph.gz) file.

    Parameters
    ----------
    input_file : str or pysam.AlignmentFile
        BAM/CRAM/fragment input.
    site_bed : str or path
        BED of sites, sorted by contig then start.  ``"-"`` reads stdin.
    chrom_sizes : str or path, optional
        ``.chrom.sizes`` file (required when ``input_file`` is a fragment file).
    output_file : str, optional
        ``.bw`` or ``.bed.gz``/`.bedGraph.gz` path.
    window_size : int, optional
        WPS window width (default 120).
    interval_size : int, optional
        Window size centered on each site midpoint (default 5000).
    min_length, max_length : int, optional
        Fragment-length filter (defaults 120/180).
    quality_threshold : int, optional
        Minimum mapping quality (default 30).
    workers : int, optional
        Worker-process count (default 1).
    verbose : bool or int, optional
        Print progress/timing.
    fraction_low, fraction_high : int, optional
        Deprecated aliases for ``min_length``/``max_length``.
    reference_file : str or Path, optional
        Reference genome (required for CRAM).

    Returns
    -------
    str or None
        The output path (or ``None`` if no output was requested).
    """
    if verbose:
        start_time = time.time()
        stderr.write(
            f"""
            Calculating aggregate WPS
            input_file: {input_file}
            site_bed: {site_bed}
            output_file: {output_file}
            window_size: {window_size}
            interval_size: {interval_size}
            quality_threshold: {quality_threshold}
            workers: {workers}
            verbose: {verbose}

            """
        )

    if input_file == "-" and site_bed == "-":
        raise ValueError("input_file and site_bed cannot both read from stdin")

    # Resolve deprecated aliases (both spellings together is an error).
    if fraction_low is not None and min_length is None:
        min_length = fraction_low
        warnings.warn(
            "fraction_low is deprecated. Use min_length instead.",
            category=DeprecationWarning,
            stacklevel=2,
        )
    elif fraction_low is not None and min_length is not None:
        warnings.warn(
            "fraction_low is deprecated. Use min_length instead.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        raise ValueError("fraction_low and min_length cannot both be specified")

    if fraction_high is not None and max_length is None:
        max_length = fraction_high
        warnings.warn(
            "fraction_high is deprecated. Use max_length instead.",
            category=DeprecationWarning,
            stacklevel=2,
        )
    elif fraction_high is not None and max_length is not None:
        warnings.warn(
            "fraction_high is deprecated. Use max_length instead.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        raise ValueError("fraction_high and max_length cannot both be specified")

    header = _read_header(input_file, chrom_sizes)
    references = [chrom for (chrom, _) in header]
    chrom_sizes_dict = dict(header)

    if verbose > 1:
        stderr.write(f"header is {header}\n")
    if verbose:
        stderr.write("Reading intervals from bed\n")

    contigs, starts, stops = _read_sites(
        site_bed, interval_size, references, chrom_sizes_dict
    )

    # BigWig requires entries in header (contig) order. BED files sorted
    # alphabetically (chr1, chr10, ... chr2) would otherwise be silently
    # dropped by pyBigWig as out-of-order writes.
    if header and contigs:
        chrom_order = {chrom: idx for idx, (chrom, _) in enumerate(header)}
        sort_indices = sorted(
            range(len(contigs)),
            key=lambda i: (chrom_order.get(contigs[i], len(header)), starts[i]),
        )
        contigs = [contigs[i] for i in sort_indices]
        starts = [starts[i] for i in sort_indices]
        stops = [stops[i] for i in sort_indices]

    try:
        chrom_sizes_intervals = [chrom_sizes_dict[contig] for contig in contigs]
    except KeyError as e:
        raise ValueError(
            f"Chrom {e} from {site_bed} is not present in {input_file} or "
            "chrom.sizes file if applicable). Please ensure that all files use "
            "the same reference genome and chromosome naming conventions."
        )

    count = len(contigs)

    if verbose:
        stderr.write("Zipping inputs\n")

    interval_list = zip(
        count * [input_file],
        contigs,
        starts,
        stops,
        chrom_sizes_intervals,
        count * [None],
        count * [window_size],
        count * [min_length],
        count * [max_length],
        count * [quality_threshold],
        count * [verbose - 2 if verbose > 2 else 0],
        count * [fraction_low],
        count * [fraction_high],
        count * [reference_file],
    )

    if verbose:
        stderr.write("Calculating wps...\n")

    pool = Pool(workers, maxtasksperchild=500)
    try:
        interval_scores = pool.imap(_wps_star, interval_list, chunksize=100)

        if isinstance(output_file, str):
            if verbose:
                stderr.write(f"Output file {output_file} specified. Opening...\n")
            if output_file.endswith(".bw"):
                _write_bigwig(output_file, header, interval_scores, stops)
            elif output_file.endswith(".bed.gz") or output_file.endswith(
                "bedGraph.gz"
            ):
                _write_bedgraph_gz(output_file, interval_scores)
            else:
                raise ValueError("output_file can only have suffix .bw")
        elif output_file is not None:
            raise TypeError(
                f'output_file is unsupported type "{type(input_file)}". '
                "output_file should be a string specifying the path of the "
                "file to output scores to."
            )
    finally:
        pool.close()

    if verbose:
        end_time = time.time()
        stderr.write(f"multi_wps took {end_time - start_time} s to complete\n")
    return output_file


def _read_header(input_file, chrom_sizes) -> list[tuple[str, int]]:
    """Derive ``(contig, length)`` pairs from a BAM/CRAM header or chrom.sizes."""
    if isinstance(input_file, pysam.AlignmentFile):
        return list(zip(input_file.references, input_file.lengths))
    if isinstance(input_file, (str, PathLike)) and str(input_file).endswith(
        (".sam", ".bam", ".cram")
    ):
        with pysam.AlignmentFile(str(input_file), "r") as bam:
            return list(zip(bam.references, bam.lengths))
    if chrom_sizes is None:
        raise ValueError("chrom_sizes must be specified for BED/Fragment files")
    return chrom_sizes_to_list(chrom_sizes)


def _read_sites(site_bed, interval_size, references, chrom_sizes_dict):
    """Parse site BED into centered, non-overlapping windows."""
    contigs: list[str] = []
    starts: list[int] = []
    stops: list[int] = []

    left_of_site = round(-interval_size / 2)
    right_of_site = round(interval_size / 2)
    assert right_of_site - left_of_site == interval_size

    bed = stdin if site_bed == "-" else open(site_bed)
    try:
        prev_contig = None
        prev_start = 0
        prev_stop = 0
        for line in bed:
            contents = line.split()
            contig = contents[0].strip()
            if int(contents[1]) > int(contents[2]):
                raise ValueError(
                    f"[multi_wps] {contig}:{contents[1]}-{contents[2]} is "
                    "invalid. Please be sure start coordinate occurs before "
                    f"stop for all intervals in {site_bed}."
                )
            if contig not in references:
                warnings.warn(
                    f"Skipping site {contig}:{int(contents[1])} from site_bed "
                    "(chrom not in chrom_sizes)",
                    UserWarning,
                )
                continue
            midpoint = (int(contents[1]) + int(contents[2])) // 2

            start = max(0, midpoint + int(left_of_site))
            stop = min(midpoint + int(right_of_site), chrom_sizes_dict[contig])

            # Truncate the previous interval where it overlaps this one.
            if contig == prev_contig and start < prev_stop:
                prev_stop = start

            if prev_contig is not None and prev_stop > prev_start:
                contigs.append(prev_contig)
                starts.append(prev_start)
                stops.append(prev_stop)

            prev_contig = contig
            prev_start = start
            prev_stop = stop

        if prev_stop > prev_start:
            contigs.append(prev_contig)
            starts.append(prev_start)
            stops.append(prev_stop)
    finally:
        if site_bed != "-":
            bed.close()

    return contigs, starts, stops


def _write_bigwig(output_file, header, interval_scores, stops) -> None:
    """Write per-interval WPS scores to a bigWig file."""
    with pbw.open(output_file, "w") as bigwig:
        bigwig.addHeader(header)
        for interval_score in interval_scores:
            contigs = interval_score["contig"]
            starts = interval_score["start"]
            scores = interval_score["wps"]

            if contigs.shape == (0,):
                continue
            try:
                bigwig.addEntries(
                    contigs[0],
                    starts[0],
                    values=scores.astype(np.float64),
                    step=1,
                    span=1,
                )
            except RuntimeError:
                stderr.write(f"{contigs[0]}:{starts[0]}-{stops[-1]}\n")
                stderr.write(
                    "/n invalid or out of order interval encountered. "
                    "Skipping to next.\n"
                )
                continue


def _write_bedgraph_gz(output_file, interval_scores) -> None:
    """Write per-position WPS scores to a gzip-compressed bedGraph."""
    with gzip.open(output_file, "wt") as bedgraph:
        for interval_score in interval_scores:
            contigs = interval_score["contig"]
            starts = interval_score["start"]
            scores = interval_score["wps"]
            stops = starts + 1

            lines = "".join(
                f"{contig}\t{start}\t{stop}\t{score}\n"
                for contig, start, stop, score in zip(contigs, starts, stops, scores)
            )
            bedgraph.write(lines)
