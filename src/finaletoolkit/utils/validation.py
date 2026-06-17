"""
Validation helpers for contig compatibility and interval bounds.
"""
from __future__ import annotations

from .logging import get_logger

logger = get_logger(__name__)

__all__ = ["validate_compatible_contigs", "valid_interval"]


def validate_compatible_contigs(
    reference_contigs: list[str] | dict[str, int],
    input_contigs: list[str] | dict[str, int],
    allow_subset: bool = True,
    validate_sizes: bool = False,
    throw_on_error: bool = True,
) -> bool:
    """Validate that input contigs are compatible with a reference.

    Parameters
    ----------
    reference_contigs : list[str] or dict[str, int]
        Reference contigs; a dict additionally provides lengths.
    input_contigs : list[str] or dict[str, int]
        Input-file contigs; a dict additionally provides lengths.
    allow_subset : bool, optional
        If ``True`` (default), input may be a subset of the reference; if
        ``False``, the contig sets must be identical.
    validate_sizes : bool, optional
        If ``True``, also require matching lengths (both args must be dicts).
    throw_on_error : bool, optional
        If ``True`` (default), raise on mismatch; otherwise log and return
        ``False``.

    Returns
    -------
    bool
        ``True`` if compatible.

    Raises
    ------
    ValueError
        If contigs are missing and ``throw_on_error`` is True.
    RuntimeError
        If sizes mismatch and ``throw_on_error`` is True.
    TypeError
        If ``validate_sizes`` is requested without dict inputs.
    """
    ref_names = set(
        reference_contigs.keys()
        if isinstance(reference_contigs, dict)
        else reference_contigs
    )
    input_names = set(
        input_contigs.keys() if isinstance(input_contigs, dict) else input_contigs
    )

    missing_in_ref = input_names - ref_names
    if missing_in_ref:
        msg = (
            "Input contains contigs not found in reference: "
            f"{sorted(missing_in_ref)}"
        )
        if throw_on_error:
            raise ValueError(msg)
        logger.error(msg)
        return False

    if not allow_subset:
        missing_in_input = ref_names - input_names
        if missing_in_input:
            msg = (
                "Reference contains contigs not found in input: "
                f"{sorted(missing_in_input)}"
            )
            if throw_on_error:
                raise ValueError(msg)
            logger.error(msg)
            return False

    if validate_sizes:
        if not isinstance(reference_contigs, dict) or not isinstance(
            input_contigs, dict
        ):
            msg = (
                "validate_sizes=True requires both reference_contigs and "
                "input_contigs to be dictionaries with lengths."
            )
            if throw_on_error:
                raise TypeError(msg)
            logger.error(msg)
            return False

        for contig in input_names:
            if reference_contigs[contig] != input_contigs[contig]:
                msg = (
                    f"Contig length mismatch for '{contig}': "
                    f"reference={reference_contigs[contig]}, "
                    f"input={input_contigs[contig]}"
                )
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
    """Validate an interval against reference contigs (and optional lengths).

    Parameters
    ----------
    reference_contigs : list[str] or dict[str, int]
        Reference contigs and optional lengths.
    contig : str
        Contig to check.
    start : int, optional
        0-based start.
    stop : int, optional
        1-based stop.
    throw_on_error : bool, optional
        If ``True``, raise on an invalid interval; otherwise log and return
        ``False`` (default ``False``).

    Returns
    -------
    bool
        ``True`` if the interval is valid.
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
            msg = (
                f"Start position {start} is out of bounds for contig "
                f"'{contig}' (length {length})."
            )
            if throw_on_error:
                raise IndexError(msg)
            logger.error(msg)
            return False
        if stop is not None and (stop < 0 or stop > length):
            msg = (
                f"Stop position {stop} is out of bounds for contig "
                f"'{contig}' (length {length})."
            )
            if throw_on_error:
                raise IndexError(msg)
            logger.error(msg)
            return False
        if start is not None and stop is not None and start >= stop:
            msg = (
                f"Invalid interval: start ({start}) must be less than stop "
                f"({stop})."
            )
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

    return True
