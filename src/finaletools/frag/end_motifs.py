from __future__ import annotations
from typing import Union

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