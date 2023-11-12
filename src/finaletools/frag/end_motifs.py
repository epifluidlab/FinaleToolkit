from __future__ import annotations
from collections.abc import Iterator
from typing import Union, Iterable
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

from finaletools.utils.utils import frag_generator
import finaletools.frag as pkg_data

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
    refseq_file: str
        A 2bit file containing the reference sequence that cfDNA
        fragments were aligned to
    quality_threshold: int, optional
        Minimum mapping quality used. Default is 30.

    """

    def __init__(
        self,
        kmer_frequencies: Iterable[tuple[str, float]],
        k: int,
        quality_threshold: int = MIN_QUALITY,
    ):
        self._dict = dict(kmer_frequencies)
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
        return self._dict.__len__()

    def __str__(self) -> str:
        return ''.join(f'{kmer}: {freq}\n' for kmer, freq in self)

    def kmers(self) -> list:
        return list(self._dict.keys())

    def frequencies(self) -> list:
        return list(self._dict.values())

    def freq(self, kmer: str) -> float:
        return self._dict[kmer]

    def to_tsv(self, output_file: str, sep: str='\t'):
        """Prints k-mer frequencies to a tsv"""
        if type(output_file) == str:
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
            raise TypeError(f'output_file must be a string.')

    def motif_diversity_score(self) -> float:
        """
        Calculates a motif diversity score (MDS) using normalized
        Shannon entropy as described by Jiang et al (2020). This
        function is generalized for any k instead of just 4-mers.
        """
        num_kmers = 4**self.k
        freq = np.array(self.frequencies())
        mds = np.sum(-freq*np.log(freq)/np.log(num_kmers))

        return mds

    @classmethod
    def from_file(
            cls,
            file_path: str,
            quality_threshold: int,
            sep: str='\t',
            header: int=0,) -> EndMotifFreqs:
        """
        Reads kmer frequency from a two-column tab-delimited file

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
            if file_path.endswith('gz'):
                is_file = True
                file = gzip.open(file_path)
            elif file_path == '-':
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
    refseq_file: str,
    k: int = 4,
    fraction_low: int = 10,
    fraction_high: int = 600,
    both_strands: bool = False,
    output_file: Union(None, str) = None,
    quality_threshold: int = MIN_QUALITY,
    verbose: Union(bool, int) = False,
) -> dict:
    """
    Function that reads fragments in the specified region from a BAM,
    SAM, or tabix indexed file and returns the 5' k-mer (default is
    4-mer) end motif counts as a structured array. This function
    reproduces the methodology of Zhou et al (2023).

    Parameters
    ----------
    input_file : str
    contig : str
    start : int
    stop : int
    refseq_file : str
    k : int, optional
    output_file : None or str, optional
    quality_threshold : int, optional
    verbose : bool or int, optional

    Return
    ------
    end_motif_freq : dict
    """
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
        refseq = py2bit.open(refseq_file, 'r')
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
                            end_motif_counts[_reverse_complement(reverse_kmer)] += 1
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


def _region_end_motifs_star(args) -> NDArray:
    results_dict = region_end_motifs(*args)
    return np.array(list(results_dict.values()), dtype='<f8')


def end_motifs(
    input_file: str,
    refseq_file: str,
    k: int = 4,
    fraction_low: int = 10,
    fraction_high: int = 600,
    both_strands: bool = False,
    output_file: Union(None, str) = None,
    quality_threshold: int = 30,
    workers: int = 1,
    verbose: Union(bool, int) = False,
) -> EndMotifFreqs:
    """
    Function that reads fragments from a BAM, SAM, or tabix indexed
    file and returns the 5' k-mer (default is 4-mer) end motif
    frequencies as a dictionary. Optionally writes data to a tsv. This
    function reproduces the methodology of Zhou et al (2023).

    Parameters
    ----------
    input_file : str
    refseq_file : str
    k : int, optional
    output_file : None or str, optional
    quality_threshold : int, optional
    workers : int, optional
    verbose : bool or int, optional

    Return
    ------
    end_motif_freq : list
    """
    if verbose:
        start_time = time()

    bases='ACGT'
    kmer_list = _gen_kmers(k, bases)

    # read chromosomes from py2bit
    try:
        refseq = py2bit.open(refseq_file, 'r')
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
