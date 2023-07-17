from __future__ import annotations
from typing import Union


def region_end_motifs(
    input_file: str,
    contig: str,
    start: int,
    stop: int,
    refseq: str,
    k: int = 4,
    output_file: Union(None, str) = None,
    verbose: Union(bool, int) = False,
) -> dict:
    """
    Function that reads fragments int the specified region from a BAM,
    SAM, or tabix indexed file and returns the 5' k-mer (default is
    4-mer) end motif frequencies as a dictionary.
    """
    return None

# XXX: implement
def end_motifs(
    input_file: str,
    refseq: str,
    k: int = 4,
    output_file: Union(None, str) = None,
    workers: int = 1,
    verbose: Union(bool, int) = False,
) -> dict:
    """
    Function that reads fragments from a BAM, SAM, or tabix indexed
    file and returns the 5' k-mer (default is 4-mer) end motif
    frequencies as a dictionary. Optionally writes data to a tsv.
    """
    return None