

from typing import Any
from .logging import get_logger

logger = get_logger(__name__)

def validate_compatible_contigs(
    reference_contigs: list[str] | dict[str, int],
    input_contigs: list[str] | dict[str, int],
    allow_subset: bool = True,
    validate_sizes: bool = False,
    throw_on_error: bool = True,
) -> bool:
    """
    Validates that the contigs in the input are compatible with the reference.

    Parameters
    ----------
    reference_contigs : list[str] | dict[str, int]
        The contigs from the reference genome. If a dict, keys are contig names
        and values are their lengths.
    input_contigs : list[str] | dict[str, int]
        The contigs from the input file (e.g., BAM, BED).
    allow_subset : bool, optional
        If True, the input contigs can be a subset of the reference contigs.
        If False, the sets must be identical. Default is True.
    validate_sizes : bool, optional
        If True, also validates that the lengths of the contigs match.
        Requires both inputs to be dictionaries. Default is False.
    throw_on_error : bool, optional
        If True, raises a ValueError or RuntimeError on mismatch.
        If False, logs the error and returns False. Default is True.

    Returns
    -------
    bool
        True if compatible, False otherwise (if throw_on_error is False).
    """
    ref_names = set(reference_contigs.keys() if isinstance(reference_contigs, dict) else reference_contigs)
    input_names = set(input_contigs.keys() if isinstance(input_contigs, dict) else input_contigs)

    # Check for missing contigs
    missing_in_ref = input_names - ref_names
    if missing_in_ref:
        msg = f"Input contains contigs not found in reference: {sorted(list(missing_in_ref))}"
        if throw_on_error:
            raise ValueError(msg)
        logger.error(msg)
        return False

    if not allow_subset:
        missing_in_input = ref_names - input_names
        if missing_in_input:
            msg = f"Reference contains contigs not found in input: {sorted(list(missing_in_input))}"
            if throw_on_error:
                raise ValueError(msg)
            logger.error(msg)
            return False

    # Validate sizes if requested
    if validate_sizes:
        if not isinstance(reference_contigs, dict) or not isinstance(input_contigs, dict):
            msg = "validate_sizes=True requires both reference_contigs and input_contigs to be dictionaries with lengths."
            if throw_on_error:
                raise TypeError(msg)
            logger.error(msg)
            return False

        for contig in input_names:
            if reference_contigs[contig] != input_contigs[contig]:
                msg = (f"Contig length mismatch for '{contig}': "
                       f"reference={reference_contigs[contig]}, input={input_contigs[contig]}")
                if throw_on_error:
                    raise RuntimeError(msg)
                logger.error(msg)
                return False

    return True

def valid_interval(
    reference_contigs: list[str] | dict[str, int],
    contig: str,
    start: int | None = None,
    stop: int | None = None,
    throw_on_error: bool = False,
) -> bool:
    """
    Validates if a given interval is valid for the reference contigs.

    Parameters
    ----------
    reference_contigs : list[str] | dict[str, int]
        The reference contigs and their optional lengths.
    contig : str
        The contig name to check.
    start : int, optional
        The start position (0-indexed).
    stop : int, optional
        The stop position (1-indexed).
    throw_on_error : bool, optional
        If True, raises a ValueError or IndexError on mismatch.
        If False, logs the error and returns False. Default is False.

    Returns
    -------
    bool
        True if the interval is valid, False otherwise.
    """
    if contig not in reference_contigs:
        msg = f"Contig '{contig}' not found in reference."
        if throw_on_error:
            raise ValueError(msg)
        logger.error(msg)
        return False

    if isinstance(reference_contigs, dict):
        length = reference_contigs[contig]
        if start is not None and (start < 0 or start >= length):
            msg = f"Start position {start} is out of bounds for contig '{contig}' (length {length})."
            if throw_on_error:
                raise IndexError(msg)
            logger.error(msg)
            return False
        if stop is not None and (stop < 0 or stop > length):
            msg = f"Stop position {stop} is out of bounds for contig '{contig}' (length {length})."
            if throw_on_error:
                raise IndexError(msg)
            logger.error(msg)
            return False
        if start is not None and stop is not None and start >= stop:
            msg = f"Invalid interval: start ({start}) must be less than stop ({stop})."
            if throw_on_error:
                raise ValueError(msg)
            logger.error(msg)
            return False
    else:
        if start is not None and start < 0:
            msg = f"Start position {start} cannot be negative."
            if throw_on_error:
                raise IndexError(msg)
            logger.error(msg)
            return False
        # Cannot check upper bounds without lengths

    return True
