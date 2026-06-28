"""
Shared infrastructure for end-motif and breakpoint-motif features.

End motifs and breakpoint motifs only differ in *where* the k-mer is read
relative to each fragment end; the container classes, motif-diversity score,
file I/O, deprecated-alias handling, and multiprocessing drivers are identical.
This module holds that shared code so the two feature modules stay thin.
"""
from __future__ import annotations

import gzip
import warnings
from collections.abc import Iterator
from importlib.resources import files
from multiprocessing import Pool
from pathlib import Path
from sys import stdin, stdout
from typing import Iterable

import numpy as np
from tqdm import tqdm

import finaletoolkit.frag as pkg_data
from finaletoolkit.utils import gen_kmers

# Path to the Zhou et al. (2023) end-motif f-profiles bundled with the package.
FPROFILE_PATH = files(pkg_data) / "data" / "end_motif_f_profiles.tsv"

# Quality threshold used by Jiang et al. (2020).
MIN_QUALITY: int = 20

_BASES = "ACGT"
_WINDOW_SIZE = 1_000_000

__all__ = ["FPROFILE_PATH", "MIN_QUALITY"]


def _normalized_shannon_mds(counts: np.ndarray, k: int) -> float:
    """Normalized-Shannon-entropy motif diversity score (Jiang et al., 2020).

    Generalized to any ``k``.  Zero-frequency motifs contribute nothing.
    """
    num_kmers = 4**k
    freq = np.asarray(counts, dtype=np.float64)
    return float(
        -np.sum(
            freq
            * np.log(
                freq,
                out=np.zeros_like(freq, dtype=np.float64),
                where=(freq != 0),
            )
            / np.log(num_kmers)
        )
    )


def resolve_motif_aliases(
    min_length: int | None,
    max_length: int | None,
    fraction_low: int | None,
    fraction_high: int | None,
) -> tuple[int | None, int | None]:
    """Resolve deprecated ``fraction_low``/``fraction_high`` aliases.

    Preserves the original behavior: using a deprecated alias warns, and
    supplying both the alias and its modern counterpart raises ``ValueError``.
    """
    if fraction_low is not None and min_length is None:
        min_length = fraction_low
        warnings.warn(
            "fraction_low is deprecated. Use min_length instead.",
            category=DeprecationWarning,
            stacklevel=3,
        )
    elif fraction_low is not None and min_length is not None:
        warnings.warn(
            "fraction_low is deprecated. Use min_length instead.",
            category=DeprecationWarning,
            stacklevel=3,
        )
        raise ValueError("fraction_low and min_length cannot both be specified")

    if fraction_high is not None and max_length is None:
        max_length = fraction_high
        warnings.warn(
            "fraction_high is deprecated. Use max_length instead.",
            category=DeprecationWarning,
            stacklevel=3,
        )
    elif fraction_high is not None and max_length is not None:
        warnings.warn(
            "fraction_high is deprecated. Use max_length instead.",
            category=DeprecationWarning,
            stacklevel=3,
        )
        raise ValueError("fraction_high and max_length cannot both be specified")

    return min_length, max_length


class _MotifFreqs:
    """Base class storing genome-wide motif k-mer frequencies.

    Subclassed by :class:`~finaletoolkit.frag.EndMotifFreqs` and
    :class:`~finaletoolkit.frag.BreakpointMotifFreqs` (which differ only in the
    default ``quality_threshold``).
    """

    def __init__(
        self,
        kmer_frequencies: Iterable[tuple[str, float]],
        k: int,
        quality_threshold: int = MIN_QUALITY,
    ) -> None:
        self.freq_dict = dict(kmer_frequencies)
        self.k = k
        self.quality_threshold = quality_threshold
        if not all(len(kmer) == k for kmer in self.freq_dict):
            raise ValueError(
                "kmer_frequencies contains a kmer with length not equal to k."
            )

    def __iter__(self) -> Iterator:
        return ((kmer, frequency) for kmer, frequency in self.freq_dict.items())

    def __len__(self) -> int:
        return len(self.freq_dict)

    def __str__(self) -> str:
        return "".join(f"{kmer}: {freq}\n" for kmer, freq in self)

    def kmers(self) -> list:
        """Return the list of k-mers."""
        return list(self.freq_dict.keys())

    def frequencies(self) -> list:
        """Return the list of frequencies."""
        return list(self.freq_dict.values())

    def freq(self, kmer: str) -> float:
        """Return the frequency of a single k-mer."""
        return self.freq_dict[kmer]

    def to_tsv(self, output_file: str | Path, sep: str = "\t") -> None:
        """Write k-mer frequencies to a TSV (``"-"`` writes to stdout)."""
        if not isinstance(output_file, (str, Path)):
            raise TypeError("output_file must be a string or path.")
        output_is_file = False
        try:
            if str(output_file) == "-":
                output = stdout
            else:
                output_is_file = True
                output = open(output_file, "w")
            for kmer, freq in self:
                output.write(f"{kmer}{sep}{freq}\n")
        finally:
            if output_is_file:
                output.close()

    def motif_diversity_score(self) -> float:
        """Normalized-Shannon-entropy motif diversity score (any ``k``)."""
        return _normalized_shannon_mds(np.array(self.frequencies()), self.k)

    @classmethod
    def from_file(
        cls,
        file_path: str | Path,
        quality_threshold: int,
        sep: str = "\t",
        header: int = 0,
    ):
        """Read k-mer frequencies from a two-column file.

        Parameters
        ----------
        file_path : str or Path
            Input file (``"-"`` reads stdin, ``.gz`` is decompressed).
        quality_threshold : int
            Stored on the resulting object.
        sep : str, optional
            Field delimiter (default tab).
        header : int, optional
            Number of leading lines to skip.
        """
        file = None
        is_file = False
        try:
            if str(file_path).endswith("gz"):
                is_file = True
                file = gzip.open(file_path, "rt")
            elif str(file_path) == "-":
                file = stdin
            else:
                is_file = True
                file = open(file_path, "rt")

            for _ in range(header):
                file.readline()

            freq_list = []
            lines = file.readlines()
            k = len(lines[header].split(sep)[0])  # infer k from first entry

            for line in lines:
                line_data = line.split(sep)
                if len(line_data) != 2:
                    break
                freq_list.append((line_data[0], float(line_data[1])))
                if k != len(line_data[0]):
                    raise RuntimeError(
                        "File contains k-mers of inconsistent length."
                    )
            if (length := len(freq_list)) != 4**k:
                raise RuntimeError(
                    f"File contains {length} {k}-mers instead of the expected "
                    f"{4**k} {k}-mers."
                )
        finally:
            if is_file and file is not None:
                file.close()
        return cls(freq_list, k, quality_threshold)


class _MotifsIntervals:
    """Base class storing motif k-mer counts stratified over intervals.

    Subclassed by :class:`~finaletoolkit.frag.EndMotifsIntervals` and
    :class:`~finaletoolkit.frag.BreakpointMotifsIntervals`.
    """

    def __init__(
        self,
        intervals: list[tuple[tuple, dict]],
        k: int,
        quality_threshold: int = MIN_QUALITY,
    ) -> None:
        self.intervals = intervals
        self.k = k
        self.quality_threshold = quality_threshold
        if not all(len(freqs) == 4**k for _, freqs in intervals):
            raise ValueError(
                "bins contains results for kmer with length not equal to k."
            )

    def __iter__(self) -> Iterator:
        return (interval for interval in self.intervals)

    def __len__(self) -> int:
        return len(self.intervals)

    def __str__(self) -> str:
        return f"{type(self).__name__} over {len(self.intervals)} intervals."

    @classmethod
    def from_file(
        cls,
        file_path: str,
        quality_threshold: int,
        sep: str = ",",
        header: int = 0,
    ):
        """Read interval-stratified motif counts from a table.

        Columns are ``contig, start, stop, name, count, <kmers...>``.
        """
        file = None
        is_file = False
        try:
            if file_path.endswith("gz"):
                is_file = True
                # Read mode (the original opened gz files in write mode here,
                # which truncated them; corrected to "rt").
                file = gzip.open(file_path, "rt")
            elif file_path == "-":
                file = stdin
            else:
                is_file = True
                file = open(file_path)

            for _ in range(header):
                file.readline()

            intervals = []
            lines = file.readlines()
            _, _, _, _, _, *kmers = lines[0].split(sep)
            k = round(np.log(len(kmers)) / np.log(4))
            assert 4**k == len(kmers), f"k={k} but should be {len(kmers)}."

            for line in lines[1:]:
                contig, start, stop, name, count, *freqs = line.split(sep)
                start, stop = int(start), int(stop)
                float_freqs = [float(freq) for freq in freqs]
                dict_freqs = dict(zip(kmers, float_freqs))
                intervals.append(((contig, start, stop, name), dict_freqs))
        finally:
            if is_file and file is not None:
                file.close()
        return cls(intervals, k, quality_threshold)

    def freq(self, kmer: str) -> dict[tuple[str, int, int], float]:
        """Map each interval to its frequency for a single k-mer."""
        return dict(
            (*interval, freq[kmer]) for interval, freq in self.intervals
        )

    def motif_diversity_score(self) -> list[tuple[tuple, float]]:
        """Regional motif diversity score (rMDS) for each region (any ``k``)."""
        mds = []
        for interval, kmers in self.intervals:
            counts = np.array(list(kmers.values()))
            total = np.sum(counts)
            try:
                region_mds = _normalized_shannon_mds(counts / total, self.k)
            except RuntimeWarning:
                region_mds = np.nan
            mds.append((interval, region_mds))
        return mds

    def mds_bed(self, output_file: str | Path, sep: str = "\t") -> None:
        """Write the regional MDS (rMDS) of each region to a BED/bedgraph file."""
        mds = self.motif_diversity_score()
        with open(output_file, "w") as out:
            for interval, region_mds in mds:
                contig, start, stop, name = interval
                row = sep.join(
                    [contig, str(start), str(stop), name, str(region_mds)]
                )
                out.write(f"{row}\n")

    def to_tsv(
        self,
        output_file: str | Path,
        calc_freq: bool = True,
        sep: str = "\t",
    ) -> None:
        """Write all intervals and motif frequencies/counts to a table.

        Columns: ``contig, start, stop, name, count, <kmers...>``.  When
        ``calc_freq`` is ``True`` each motif value is a 6-decimal frequency
        (``NaN`` for zero-count intervals); otherwise raw counts are written.
        """
        if not isinstance(output_file, (str, Path)):
            raise TypeError("output_file must be a string or path.")
        output_is_file = False
        try:
            if str(output_file) == "-":
                output = stdout
            else:
                output_is_file = True
                output = open(output_file, "w")

            kmers = gen_kmers(self.k, _BASES)
            output.write(
                sep.join(["contig", "start", "stop", "name", "count", *kmers])
            )
            output.write("\n")
            for interval, freqs in self.intervals:
                count = sum(freqs.values())
                if calc_freq:
                    values = [
                        f"{(freq / count):.6f}" if count != 0 else "NaN"
                        for freq in freqs.values()
                    ]
                else:
                    values = [str(freq) for freq in freqs.values()]
                output.write(
                    sep.join(
                        [
                            interval[0],
                            str(interval[1]),
                            str(interval[2]),
                            str(interval[3]),
                            str(count),
                            *values,
                        ]
                    )
                )
                output.write("\n")
        finally:
            if output_is_file:
                output.close()

    def _to_record(
        self, kmer: str, output_file, calc_freq: bool, sep: str, include_name: bool
    ) -> None:
        """Shared writer for :meth:`to_bed`/:meth:`to_bedgraph`."""
        if not isinstance(output_file, (str, Path)):
            raise TypeError("output_file must be a string.")
        output_is_file = False
        try:
            if str(output_file) == "-":
                output = stdout
            else:
                output_is_file = True
                output = open(output_file, "w")

            for interval, freqs in self.intervals:
                count = sum(freqs.values())
                if calc_freq:
                    value = f"{(freqs[kmer] / count):.6f}" if count != 0 else "NaN"
                else:
                    value = freqs[kmer]
                fields = [interval[0], str(interval[1]), str(interval[2])]
                if include_name:
                    fields.append(interval[3])
                fields.append(value)
                output.write(sep.join(fields))
                output.write("\n")
        finally:
            if output_is_file:
                output.close()

    def to_bedgraph(
        self,
        kmer: str,
        output_file: str | Path,
        calc_freq: bool = True,
        sep: str = "\t",
    ) -> None:
        """Write one k-mer's per-interval value to a bedGraph."""
        self._to_record(kmer, output_file, calc_freq, sep, include_name=False)

    def to_bed(
        self,
        kmer: str,
        output_file: str | Path,
        calc_freq: bool = True,
        sep: str = "\t",
    ) -> None:
        """Write one k-mer's per-interval value to a BED."""
        self._to_record(kmer, output_file, calc_freq, sep, include_name=True)


# -- multiprocessing drivers shared by both feature types -------------------


def _genome_window_args(
    input_file,
    refseq_file,
    chroms,
    k,
    min_length,
    max_length,
    both_strands,
    negative_strand,
    quality_threshold,
    verbose,
) -> list[tuple]:
    """Build per-1Mb-window argument tuples covering every contig."""
    intervals = []
    for chrom, chrom_length in chroms.items():
        for start in range(0, chrom_length - _WINDOW_SIZE, _WINDOW_SIZE):
            intervals.append(
                (
                    input_file,
                    chrom,
                    start,
                    start + _WINDOW_SIZE,
                    refseq_file,
                    k,
                    min_length,
                    max_length,
                    both_strands,
                    negative_strand,
                    None,
                    quality_threshold,
                    verbose - 2 if verbose > 2 else 0,
                )
            )
        intervals.append(
            (
                input_file,
                chrom,
                chrom_length - chrom_length % _WINDOW_SIZE,
                chrom_length,
                refseq_file,
                k,
                min_length,
                max_length,
                both_strands,
                negative_strand,
                None,
                quality_threshold,
                verbose - 2 if verbose > 2 else 0,
            )
        )
    return intervals


def aggregate_genome_motifs(
    intervals,
    region_star,
    freqs_class,
    k,
    quality_threshold,
    workers,
    verbose,
    read_desc,
    count_desc,
):
    """Fan window args out to ``region_star`` and sum into a freqs object."""
    pool = Pool(workers)
    try:
        counts_iter = pool.imap(
            region_star,
            tqdm(intervals, read_desc, position=0) if verbose else intervals,
            chunksize=min(int(len(intervals) / workers / 2 + 1), 1000),
        )
        ccounts = np.zeros((4**k,), np.float64)
        for count in (
            tqdm(counts_iter, count_desc, len(intervals), position=1)
            if verbose
            else counts_iter
        ):
            ccounts = ccounts + count
    finally:
        pool.close()

    frequencies = ccounts / np.sum(ccounts)
    return freqs_class(zip(gen_kmers(k, _BASES), frequencies), k, quality_threshold)


def parse_intervals_arg(intervals) -> list[tuple]:
    """Coerce a BED path or list-of-tuples into ``(chrom, start, stop, name)``."""
    if type(intervals) is str:
        with open(intervals, "r") as interval_file:
            return [
                (
                    chrom,
                    int(start),
                    int(stop),
                    name[0] if len(name) > 0 else ".",
                )
                for chrom, start, stop, *name in (
                    line.split() for line in interval_file.readlines()
                )
            ]
    if isinstance(intervals, list):
        return intervals
    raise TypeError("Intervals should be string or list.")


def aggregate_interval_motifs(
    input_file,
    refseq_file,
    intervals,
    region_dict_star,
    intervals_class,
    k,
    min_length,
    max_length,
    both_strands,
    negative_strand,
    quality_threshold,
    workers,
    verbose,
    read_desc,
):
    """Fan interval args out to ``region_dict_star`` and build a class instance."""
    intervals_tuples = parse_intervals_arg(intervals)
    mp_intervals = [
        (
            input_file,
            chrom,
            start,
            stop,
            refseq_file,
            k,
            min_length,
            max_length,
            both_strands,
            negative_strand,
            None,
            quality_threshold,
            verbose - 2 if verbose > 2 else 0,
        )
        for chrom, start, stop, *_ in intervals_tuples
    ]

    pool = Pool(workers)
    try:
        counts_iter = pool.imap(
            region_dict_star,
            tqdm(mp_intervals, read_desc, position=0) if verbose else mp_intervals,
            chunksize=min(int(len(mp_intervals) / workers / 2 + 1), 1000),
        )
        result = intervals_class(
            [
                (interval, counts)
                for interval, counts in zip(intervals_tuples, counts_iter)
            ],
            k,
            quality_threshold,
        )
    finally:
        pool.close()
    return result


def write_motif_freqs(results, output_file: str) -> None:
    """Write a freqs object to TSV/CSV based on the output suffix."""
    if output_file is None:
        return
    if output_file.endswith(".csv"):
        results.to_tsv(output_file, sep=",")
    else:
        results.to_tsv(output_file)
