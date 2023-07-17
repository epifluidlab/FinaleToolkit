from __future__ import annotations
from typing import Union

import py2bit
import numpy as np
from numpy.typing import NDArray

from finaletools.utils.utils import frag_array

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
    output_file: Union(None, str) = None,
    quality_threshold: int = 30,
    verbose: Union(bool, int) = False,
) -> dict:
    """
    Function that reads fragments int the specified region from a BAM,
    SAM, or tabix indexed file and returns the 5' k-mer (default is
    4-mer) end motif counts as a structured array.

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
    # numpy array of fragments
    frag_ends = frag_array(
        input_file,
        contig,
        quality_threshold,
        start,
        stop,
        fraction_low=4,
        fraction_high=100000,
    )
    # create dict where keys are kmers and values are counts
    bases='ACGT'
    kmer_list = _gen_kmers(k, bases)
    end_motif_counts = dict(zip(kmer_list, 4**k*[0]))

    # TODO: accept other reference file types e.g. FASTA
    # count end motifs
    try:
        refseq = py2bit.open(refseq_file, 'r')
        for frag in frag_ends:
            # forward end-motif
            forward_kmer = refseq.sequence(
                contig, int(frag[0]), int(frag[0]+k)
            )
            if 'N' not in forward_kmer:
                end_motif_counts[forward_kmer] += 1

            # reverse end-motif
            reverse_kmer = refseq.sequence(
                contig, int(frag[1]-k), int(frag[1])
            )
            if 'N' not in reverse_kmer:
                end_motif_counts[_reverse_complement(reverse_kmer)] += 1
    finally:
        refseq.close()

    return end_motif_counts

# XXX: implement
def end_motifs(
    input_file: str,
    refseq_file: str,
    k: int = 4,
    output_file: Union(None, str) = None,
    quality_threshold: int = 30,
    workers: int = 1,
    verbose: Union(bool, int) = False,
) -> NDArray:
    """
    Function that reads fragments from a BAM, SAM, or tabix indexed
    file and returns the 5' k-mer (default is 4-mer) end motif
    frequencies as a dictionary. Optionally writes data to a tsv.

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
    return None