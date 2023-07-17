from __future__ import annotations
from typing import Union

from finaletools.utils.utils import frag_array

# XXX: implement
def region_end_motifs(
    input_file: str,
    contig: str,
    start: int,
    stop: int,
    refseq: str,
    k: int = 4,
    output_file: Union(None, str) = None,
    quality_threshold: int = 30,
    verbose: Union(bool, int) = False,
) -> dict:
    """
    Function that reads fragments int the specified region from a BAM,
    SAM, or tabix indexed file and returns the 5' k-mer (default is
    4-mer) end motif frequencies as a dictionary.

    Parameters
    ----------
    input_file : str
    contig : str
    start : int
    stop : int
    refseq : str
    k : int, optional
    output_file : None or str, optional
    quality_threshold : int, optional
    verbose : bool or int, optional

    Return
    ------
    end_motif_freq : dict
    """
    return None

# XXX: implement
def end_motifs(
    input_file: str,
    refseq: str,
    k: int = 4,
    output_file: Union(None, str) = None,
    quality_threshold: int = 30,
    workers: int = 1,
    verbose: Union(bool, int) = False,
) -> dict:
    """
    Function that reads fragments from a BAM, SAM, or tabix indexed
    file and returns the 5' k-mer (default is 4-mer) end motif
    frequencies as a dictionary. Optionally writes data to a tsv.

    Parameters
    ----------
    input_file : str
    refseq : str
    k : int, optional
    output_file : None or str, optional
    quality_threshold : int, optional
    workers : int, optional
    verbose : bool or int, optional

    Return
    ------
    end_motif_freq : dict
    """
    return None