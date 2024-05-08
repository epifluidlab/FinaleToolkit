from __future__ import annotations
from collections.abc import Iterator
from typing import Union, Iterable, Tuple
from multiprocessing import Pool
from time import time
from sys import stderr, stdout, stdin
import gzip
try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files
from pathlib import Path
from collections import UserDict

import tqdm
import py2bit
import numpy as np
from numpy.typing import NDArray

from finaletoolkit.utils.utils import frag_generator, _get_intervals
import finaletoolkit.frag as pkg_data

# path to tsv containing f-profiles from Zhou et al (2023)
FPROFILE_PATH: Path = (files(pkg_data) / 'data' / 'end_motif_f_profiles.tsv')

# quality threshold used by Jiang et al (2020)
MIN_QUALITY: int = 20

class EndMotifFreqs():
    """
    Class that stores frequencies of end-motif k-mer frequencies and
    contains methods to manipulate this data.

    Parameters
    ----------
    kmer_frequencies : Iterable
        A Iterable of tuples, each containing a str representing a k-mer
        and a float representing its frequency
    k : int
        Size of k-mers
    quality_threshold: int, optional
        Minimum mapping quality used. Default is 30.
    """

    def __init__(
        self,
        kmer_frequencies: Iterable[tuple[str, float]],
        k: int,
        quality_threshold: int = MIN_QUALITY,
    ):
        self.freq_dict = dict(kmer_frequencies)
        self.k = k
        self.quality_threshold = quality_threshold
        if not all(len(kmer) == k for kmer, _ in kmer_frequencies):
            raise ValueError(
                'kmer_frequencies contains a kmer with length not equal'
                ' to k.'
            )

    def __iter__(self) -> Iterator:
        return ((kmer, frequency)
                for (kmer, frequency)
                in zip(self.kmers(), self.frequencies()))

    def __len__(self) -> int:
        return self.freq_dict.__len__()

    def __str__(self) -> str:
        return ''.join(f'{kmer}: {freq}\n' for kmer, freq in self)

    def kmers(self) -> list:
        return list(self.freq_dict.keys())

    def frequencies(self) -> list:
        return list(self.freq_dict.values())

    def freq(self, kmer: str) -> float:
        return self.freq_dict[kmer]

    def to_tsv(self, output_file: Union[str,Path], sep: str='\t'):
        """Prints k-mer frequencies to a tsv"""
        if isinstance(output_file, str) or isinstance(output_file, Path):
            try:
                # open file based on name
                output_is_file = False
                if output_file == '-':
                    output = stdout
                else:
                    output_is_file = True
                    output = open(output_file, 'w')

                # write to file
                for kmer, freq in self:
                    output.write(f'{kmer}{sep}{freq}\n')

            finally:
                if output_is_file:
                    output.close()
        else:
            raise TypeError(f'output_file must be a string or path.')

    def motif_diversity_score(self) -> float:
        """
        Calculates a motif diversity score (MDS) using normalized
        Shannon entropy as described by Jiang et al (2020). This
        function is generalized for any k instead of just 4-mers.
        """
        num_kmers = 4**self.k
        freq = np.array(self.frequencies())
        # if freq is 0, ignore
        mds = -np.sum(
            freq * np.log(
                freq, out=np.zeros_like(freq, dtype=np.float64), where=(freq!=0)
                ) / np.log(num_kmers))
        return mds

    @classmethod
    def from_file(
            cls,
            file_path: Union[str, Path],
            quality_threshold: int,
            sep: str='\t',
            header: int=0,) -> EndMotifFreqs:
        """
        Reads kmer frequency from a two-column tab-delimited file.

        Parameters
        ---------
        file_path : str
            Path string containing path to file.
        sep : str, optional
            Delimiter used in file.
        header : int, optional
            Number of lines to ignore at the head of the file.

        Return
        ------
        kmer_freqs : EndMotifFreqs
        """
        try:
            # open file
            is_file = False
            if str(file_path).endswith('gz'):
                is_file = True
                file = gzip.open(file_path)
            elif str(file_path) == '-':
                file = stdin
            else:
                is_file = True
                file = open(file_path)

            # ignore header
            for _ in range(header):
                file.readline()

            freq_list = []
            lines = file.readlines()
            line = lines[header].split(sep)
            k = len(line[0])    # infer k from first entry

            for line in lines:
                line_data = line.split(sep)
                if len(line_data) != 2:
                    break
                freq_list.append((line_data[0], float(line_data[1])))
                if k != len(line_data[0]):
                    raise RuntimeError(
                        'File contains k-mers of inconsistent length.'
                    )
            if length := len(freq_list) != 4**k:
                raise RuntimeError(
                    f'File contains {length} {k}-mers instead of the expected'
                    f' {4**k} {k}-mers.'
                )
        finally:
            if is_file:
                file.close()
        return cls(freq_list, k, quality_threshold)


class EndMotifsIntervals():
    """
    Class that stores frequencies of end-motif k-mers over
    user-specified intervals and contains methods to manipulate this
    data.

    Parameters
    ----------
    intervals : Iterable
        A collection of tuples, each containing a tuple representing
        a genomic interval (chrom, 0-based start, 1-based stop) and a
        dict that maps kmers to frequencies in the interval.
    k : int
        Size of k-mers
    quality_threshold: int, optional
        Minimum mapping quality used. Default is 30.
    """

    def __init__(
        self,
        intervals: Iterable[tuple[tuple, dict]],
        k: int,
        quality_threshold: int = MIN_QUALITY,
    ):
        self.intervals = intervals
        self.k = k
        self.quality_threshold = quality_threshold
        if not all(len(freqs) == 4**k for _, freqs in intervals):
            raise ValueError(
                'bins contains results for kmer with length not equal'
                ' to k.'
            )

    def __iter__(self) -> Iterator:
        return (interval for interval in self.intervals)

    def __len__(self) -> int:
        return self.intervals.__len__()

    def __str__(self) -> str:
        return f'EndMotifsIntervals over {len(self.intervals)} intervals.'
    
    @classmethod
    def from_file(
            cls,
            file_path: str,
            quality_threshold: int,
            sep: str = ',',) -> EndMotifFreqs:
        """
        Reads kmer frequency from a tab-delimited file. Expected columns
        are contig, start, stop, name, count, *kmers. Because
        exporting to file includes an option to turn counts to a fraction,
        this doesn't perfectly correspond to replicating the other file.

        Parameters
        ---------
        file_path : str
            Path string containing path to file.
        quality_threshold : int
            MAPQ filter used. Only used for some calculations.
        sep : str, optional
            Delimiter used in file.

        Return
        ------
        kmer_freqs : EndMotifFreqs
        """
        try:
            # open file
            is_file = False
            if file_path.endswith('gz'):
                is_file = True
                file = gzip.open(file_path)
            elif file_path == '-':
                file = stdin
            else:
                is_file = True
                file = open(file_path)

            intervals = []
            lines = file.readlines()
            _,_,_,_,_,*kmers = lines[0].split(sep)
            k = round(np.log(len(kmers))/np.log(4))
            assert 4**k == len(kmers), f"k={k} but should be {len(kmers)}."

            for line in lines[1:]:
                contig, start, stop, name, count, *freqs = line.split(sep)
                start, stop = int(start), int(stop)
                float_freqs = [float(freq) for freq in freqs]
                dict_freqs = dict(zip(kmers, float_freqs))
                intervals.append(((contig, start, stop, name), dict_freqs))
        finally:
            if is_file:
                file.close()
        return cls(intervals, k, quality_threshold)

    def freq(self, kmer: str) -> list[Tuple[str, int, int, float]]:
        """
        Returns a list of intervals and associated frquency for given
        kmer. Results are in the form (chrom, 0-based start, 1-based
        stop, frequency).
        """
        return dict(
            [(*interval, freq[kmer]) for interval, freq in self.intervals]
        )

    def motif_diversity_score(self) -> list[tuple[tuple, float]]:
        """
        Calculates a motif diversity score (MDS) for each interval using
        normalized Shannon entropy as described by Jiang et al (2020). This
        function is generalized for any k instead of just 4-mers.
        """
        num_kmers = 4**self.k
        mds = []
        for interval, kmers in self.intervals:
            counts = np.array(list(kmers.values()))
            freq = counts / np.sum(counts)
            try:
                interval_mds = -np.sum(
                    freq * np.log(
                        freq, out=np.zeros_like(freq, dtype=np.float64), where=(freq!=0)
                        ) / np.log(num_kmers))
            except RuntimeWarning:
                interval_mds = np.NaN
            mds.append((interval, interval_mds))
        return mds

    def mds_bed(self, output_file: Union[str, Path], sep: str='\t'):
        """Writes MDS for each interval to a bed/bedgraph file."""
        mds = self.motif_diversity_score()
        with open(output_file, 'w') as out:
            for interval, interval_mds in mds:
                contig, start, stop, name = interval
                temp_str = sep.join(
                        [contig, str(start), str(stop), name, str(interval_mds)]
                    )
                out.write(
                    f"{temp_str}\n"
                )
    
    def to_tsv(
            self,
            output_file: Union[str, Path],
            calc_freq: bool=True, sep: str='\t'):
        """
        Writes all intervals and associated frquencies to file. Columns
        are contig, start, stop, name, count, *kmers.
        
        Parameters
        ----------
        output_file: str
            File to write frequencies to.
        calc_freq: bool, optional
            Calculates frequency of motifs if true. Otherwise, writes counts
            for each motif. Default is true.
        sep: str, optional
            Separator for table. Tab-separated by default.
        """
        if isinstance(output_file, str) or isinstance(output_file, Path):
            try:
                # open file based on name
                output_is_file = False
                if str(output_file) == '-':
                    output = stdout
                else:
                    output_is_file = True
                    output = open(output_file, 'w')

                # write to file
                kmers = _gen_kmers(self.k, 'ACGT')
                # header
                output.write(sep.join(['contig','start','stop','name','count',*kmers]))
                output.write('\n')
                # data
                for interval, freqs in self.intervals:
                    count = sum(freqs.values())
                    output.write(
                        sep.join([
                            interval[0],
                            str(interval[1]),
                            str(interval[2]),
                            str(interval[3]),
                            str(count), 
                            *([str(freq) for freq in freqs.values()]
                             if not calc_freq
                             else [f"{(freq/count):.6f}"
                                   if count!=0
                                   else "NaN" 
                                for freq
                                in freqs.values()])
                        ])
                    )
                    output.write('\n')

            finally:
                if output_is_file:
                    output.close()
        else:
            raise TypeError(f'output_file must be a string or path.')
    
    def to_bedgraph(
            self,
            kmer: str,
            output_file: Union[str, Path],
            calc_freq: bool=True,
            sep: str='\t'
        ):
        """
        Take frequency of specified kmer and writes to bedgraph.
        
        Parameters
        ----------
        output_file: str
            File to write frequencies to.
        calc_freq: bool, optional
            Calculates frequency of motifs if true. Otherwise, writes counts
            for each motif. Default is true.
        sep: str, optional
            Separator for table. Tab-separated by default.
        """
        if isinstance(output_file, str) or isinstance(output_file, Path):
            try:
                # open file based on name
                output_is_file = False
                if str(output_file) == '-':
                    output = stdout
                else:
                    output_is_file = True
                    output = open(output_file, 'w')

                # write to file
                for interval, freqs in self.intervals:
                    count = sum(freqs.values())
                    output.write(
                        sep.join([
                            interval[0],
                            str(interval[1]),
                            str(interval[2]),
                            (self.freq[kmer] if not calc_freq
                             else f"{(freqs[kmer]/count):.6f}"
                                if count!=0
                                else "NaN"
                            )
                        ])
                    )
                    output.write('\n')

            finally:
                if output_is_file:
                    output.close()
        else:
            raise TypeError(f'output_file must be a string.')
        
    def to_bed(
            self,
            kmer: str,
            output_file: Union[str, Path],
            calc_freq: bool=True,
            sep: str='\t'
        ):
        """
        Take frequency of specified kmer and writes to BED.
        
        Parameters
        ----------
        output_file: str
            File to write frequencies to.
        calc_freq: bool, optional
            Calculates frequency of motifs if true. Otherwise, writes counts
            for each motif. Default is true.
        sep: str, optional
            Separator for table. Tab-separated by default.
        """
        if isinstance(output_file, str) or isinstance(output_file, Path):
            try:
                # open file based on name
                output_is_file = False
                if str(output_file) == '-':
                    output = stdout
                else:
                    output_is_file = True
                    output = open(output_file, 'w')

                # write to file
                for interval, freqs in self.intervals:
                    count = sum(freqs.values())
                    output.write(
                        sep.join([
                            interval[0],
                            str(interval[1]),
                            str(interval[2]),
                            interval[3],
                            (self.freq[kmer] if not calc_freq
                             else f"{(freqs[kmer]/count):.6f}"
                                if count!=0
                                else "NaN"
                            )
                        ])
                    )
                    output.write('\n')

            finally:
                if output_is_file:
                    output.close()
        else:
            raise TypeError(f'output_file must be a string.')


def _gen_kmers(k: int, bases: str) -> list:
        """Function to recursively create a list of k-mers."""
        if k == 1:
            return [base for base in bases]
        else:
            kmers = []
            for k_minus_mer in _gen_kmers(k-1, bases):
                for base in bases:
                    kmers.append(k_minus_mer+base)
            return kmers


def _reverse_complement(kmer: str) -> str:
    reversed = kmer[-1::-1]
    pair_dict = {'A': 'T',
                 'T': 'A',
                 'G': 'C',
                 'C': 'G'}
    complemented = ''.join(pair_dict[base] for base in reversed)
    return complemented


def region_end_motifs(
    input_file: str,
    contig: str,
    start: int,
    stop: int,
    refseq_file: Union[str, Path],
    k: int = 4,
    fraction_low: int = 10,
    fraction_high: int = 600,
    both_strands: bool = True,
    output_file: Union[None, str] = None,
    quality_threshold: int = MIN_QUALITY,
    verbose: Union[bool, int] = False,
) -> dict:
    """
    Function that reads fragments in the specified region from a BAM,
    SAM, or tabix indexed file and returns the 5' k-mer (default is
    4-mer) end motif counts as a dictionary. This function
    reproduces the methodology of Zhou et al (2023).

    Parameters
    ----------
    input_file : str
        Path of SAM, BAM, CRAM, or Frag.gz containing pair-end reads.
    contig : str
        Name of contig or chromosome for region.
    start : int
        0-based start coordinate.
    stop : int
        1-based end coordinate.
    refseq_file : str or Path
        2bit file with reference sequence `input_file` was aligned to.
    k : int, optional
        Length of end motif kmer. Default is 4.
    fraction_low: int, optional
        Minimum fragment length.
    fraction_high: int, optional
        Maximum fragment length.
    both_strands: bool, optional
        Choose whether to use forward 5' ends only or use 5' ends for
        both ends of PE reads.
    output_file : None or str, optional
    quality_threshold : int, optional
    verbose : bool or int, optional

    Return
    ------
    end_motif_freq : dict
    """
    # NOTE: consider renaming to interval_end_motif

    if verbose:
        start_time = time()

    # iterable of fragments
    frag_ends = frag_generator(
        input_file,
        contig,
        quality_threshold,
        start,
        stop,
        fraction_low=fraction_low,
        fraction_high=fraction_high,
    )
    # create dict where keys are kmers and values are counts
    bases='ACGT'
    kmer_list = _gen_kmers(k, bases)
    end_motif_counts = dict(zip(kmer_list, 4**k*[0]))

    # TODO: accept other reference file types e.g. FASTA
    # count end motifs
    try:
        refseq = py2bit.open(str(refseq_file), 'r')
        if both_strands:   # both strands of fragment
            for frag in frag_ends:
                # py2bit uses 0-based for start, 1-based for end
                # forward end-motif
                forward_kmer = refseq.sequence(
                    contig, int(frag[1]), int(frag[1]+k)
                )
                assert len(forward_kmer) == k    

                if 'N' not in forward_kmer:
                    end_motif_counts[forward_kmer] += 1
                    
                # reverse end-motif
                try:
                    reverse_kmer = refseq.sequence(
                        contig, int(frag[2]-k), int(frag[2])
                    )
                    assert len(reverse_kmer) == k

                    if 'N' not in reverse_kmer:
                        end_motif_counts[_reverse_complement(reverse_kmer)] += 1
                except RuntimeError:
                    if verbose > 1:
                        stderr.write(
                            f'Attempt to read interval at {contig}:'
                            f'{int(frag[2]-k)}-{int(frag[2])} failed.'
                            'Skipping.')
                    continue
        else:
            for frag in frag_ends:
                if frag[3]: # is on forward strand or not
                    # py2bit uses 0-based for start, 1-based for end
                    # forward end-motif
                    forward_kmer = refseq.sequence(
                        contig, int(frag[1]), int(frag[1]+k)
                    )
                    assert len(forward_kmer) == k    

                    if 'N' not in forward_kmer:
                        end_motif_counts[forward_kmer] += 1
                    
                else:
                    # reverse end-motif
                    try:
                        reverse_kmer = refseq.sequence(
                            contig, int(frag[2]-k), int(frag[2])
                        )
                        assert len(reverse_kmer) == k

                        if 'N' not in reverse_kmer:
                            rc_reverse_kmer = _reverse_complement(reverse_kmer)
                            end_motif_counts[rc_reverse_kmer] += 1
                    except RuntimeError:
                        if verbose > 1:
                            stderr.write(
                                f'Attempt to read interval at {contig}:'
                                f'{int(frag[2]-k)}-{int(frag[2])} failed.'
                                'Skipping.')
                        continue

    finally:
        refseq.close()

    if verbose:
        stop_time = time()
        stderr.write(
            f'region_end_motifs took {stop_time-start_time} seconds to run\n'
        )

    return end_motif_counts


def _region_end_motifs_star(args) -> dict:
    results_dict = region_end_motifs(*args)
    return np.array(list(results_dict.values()), dtype='<f8')


def _region_end_motifs_dict_star(args) -> dict:
    results_dict = region_end_motifs(*args)
    return results_dict


def end_motifs(
    input_file: str,
    refseq_file: Union[str, Path],
    k: int = 4,
    fraction_low: int = 10,
    fraction_high: int = 600,
    both_strands: bool = False,
    output_file: Union[None, str] = None,
    quality_threshold: int = 30,
    workers: int = 1,
    verbose: Union[bool, int] = False,
) -> EndMotifFreqs:
    """
    Function that reads fragments from a BAM, SAM, or tabix indexed
    file and returns the 5' k-mer (default is 4-mer) end motif
    frequencies as a dictionary. Optionally writes data to a tsv. This
    function reproduces the methodology of Zhou et al (2023).

    Parameters
    ----------
    input_file : str
        SAM, BAM, CRAM, or Frag.gz file with paired-end reads.
    refseq_file : str or Path
        2bit file with sequence of reference genome input_file is
        aligned to.
    k : int, optional
        Length of end motif kmer. Default is 4.
    output_file : None or str, optional
        File path to write results to. Either tsv or csv.
    quality_threshold : int, optional
        Minimum MAPQ to filter.
    workers : int, optional
        Number of worker processes.
    verbose : bool or int, optional

    Return
    ------
    end_motif_freq : EndMotifFreqs
    """
    if verbose:
        start_time = time()

    bases='ACGT'
    kmer_list = _gen_kmers(k, bases)

    # read chromosomes from py2bit
    try:
        refseq = py2bit.open(str(refseq_file), 'r')
        chroms: dict = refseq.chroms()
    finally:
        refseq.close()

    # generate list of inputs
    intervals = []
    window_size = 1000000
    for chrom, chrom_length in chroms.items():
        for start in range(0, chrom_length, window_size):
            intervals.append((
                input_file,
                chrom,
                start,
                start+window_size,
                refseq_file,
                k,
                fraction_low,
                fraction_high,
                both_strands,
                None,
                quality_threshold,
                verbose - 2 if verbose > 2 else 0
            ))
        intervals.append((
            input_file,
            chrom,
            chrom_length - chrom_length%window_size,
            chrom_length,
            refseq_file,
            k,
            fraction_low,
            fraction_high,
            both_strands,
            None,
            quality_threshold,
            verbose - 2 if verbose > 2 else 0
        ))

    # use process pool to count kmers
    try:
        # open pool
        pool = Pool(workers)

        # uses tqdm loading bar if verbose == True
        counts_iter = pool.imap(
            _region_end_motifs_star,
            tqdm.tqdm(intervals, 'Reading 1mb windows', position=0)if verbose else intervals,
            chunksize=min(int(len(intervals)/workers/2+1), 1000)
        )

        ccounts = np.zeros((4**k,), np.float64)
        for count in tqdm.tqdm(counts_iter, 'Counting end-motifs', len(intervals), position=1) if verbose else counts_iter:
            ccounts = ccounts + count

    finally:
        pool.close()

    frequencies = ccounts/np.sum(ccounts)

    results = EndMotifFreqs(
        zip(kmer_list, frequencies),
        k,
        quality_threshold,
    )

    if output_file is not None:
        if output_file.endswith('.csv'):
            results.to_tsv(output_file, sep=',')
        else:
            results.to_tsv(output_file)

    if verbose:
        stop_time = time()
        tqdm.tqdm.write(
            f'end_motifs took {stop_time-start_time} seconds to run\n'
        )

    return results


def interval_end_motifs(
    input_file: str,
    refseq_file: Union[str, Path],
    intervals: Union[str, Iterable[Tuple[str,int,int,str]]],
    k: int = 4,
    fraction_low: int = 10,
    fraction_high: int = 600,
    both_strands: bool = True,
    output_file: Union[None, str] = None,
    quality_threshold: int = 30,
    workers: int = 1,
    verbose: Union[bool, int] = False,
) -> EndMotifsIntervals:
    """
    Function that reads fragments from a BAM, SAM, or tabix indexed
    file and user-specified intervals and returns the 5' k-mer
    (default is 4-mer) end motif. Optionally writes data to a tsv.

    Parameters
    ----------
    input_file : str
        Path of SAM, BAM, CRAM, or Frag.gz containing pair-end reads.
    refseq_file : str or Path
        Path of 2bit file for reference genome that reads are aligned to.
    intervals : str or tuple
        Path of BED file containing intervals or list of tuples
        (chrom, start, stop, name).
    k : int, optional
        Length of end motif kmer. Default is 4.
    output_file : None or str, optional
        File path to write results to. Either tsv or csv.
    quality_threshold : int, optional
        Minimum MAPQ to filter.
    workers : int, optional
        Number of worker processes.
    verbose : bool or int, optional

    Return
    ------
    end_motif_freq : EndMotifIntervals
    """
    if verbose:
        start_time = time()

    bases='ACGT'
    kmer_list = _gen_kmers(k, bases)

    # generate list of inputs
    if type(intervals) is str:
        with open(intervals, 'r') as interval_file:
            intervals_tuples = [    # parses file lines into a tuple generator
                (chrom, int(start), int(stop),
                    name[0] if len(name) > 0 else '.')
                for chrom, start, stop, *name
                in (line.split() for line in interval_file.readlines())
            ]
    elif isinstance(intervals, Iterable):
        intervals_tuples = intervals
    else:
        raise TypeError("Intervals should be string or list.")
    
    mp_intervals = []   # args to be fed into pool processes
    for chrom, start, stop, *_ in intervals_tuples:
        mp_intervals.append((
            input_file,
            chrom,
            start,
            stop,
            refseq_file,
            k,
            fraction_low,
            fraction_high,
            both_strands,
            None,
            quality_threshold,
            verbose - 2 if verbose > 2 else 0
        )
    )

    # use process pool to aggregate kmers
    try:
        # open pool
        pool = Pool(workers)

        # uses tqdm loading bar if verbose == True
        counts_iter = pool.imap(
            _region_end_motifs_dict_star,
            tqdm.tqdm(
                mp_intervals,
                'Reading intervals',
                position=0) if verbose else mp_intervals,
            chunksize=min(int(len(intervals)/workers/2+1), 1000)
        )

    finally:
        pool.close()
    results = EndMotifsIntervals(
        [(interval, counts)
         for interval, counts
         in zip(intervals_tuples, counts_iter)],
         k,
         quality_threshold,
    )

    if output_file is not None:
        if output_file.endswith('.csv'):
            results.to_tsv(output_file, sep=',')
        else:
            results.to_tsv(output_file)

    if verbose:
        stop_time = time()
        tqdm.tqdm.write(
            f'end_motifs took {stop_time-start_time} seconds to run\n'
        )

    return results


def _cli_mds(
    file_path: str,
    sep: str = '\t',
    header: int = 0,
) -> float:
    """Function for commandline acces to MDS from a tsv file."""
    # 30 is used as a placeholder for the quality threshold. It is not
    # used to calculate MDS and can be ignored.
    motifs = EndMotifFreqs.from_file(
        file_path,
        30,
        sep,
        header,
    )
    mds = motifs.motif_diversity_score()
    stdout.write(f'{mds}\n')

def _cli_interval_mds(
    file_path: str,
    file_out: str,
    sep: str = ',',
    header: int = 0,
) -> float:
    # 30 is used as a placeholder for the quality threshold. It is not
    # used to calculate MDS and can be ignored.
    motifs = EndMotifsIntervals.from_file(
        file_path,
        30,
        sep,
    )
    motifs.mds_bed(file_out)
