"""
None-tolerant comparison operators.

These treat a ``None`` operand as "unbounded", returning ``True``.  They make
optional length/region bounds easy to express without branching:
``_none_geq(length, min_length)`` is ``True`` when ``min_length`` is ``None``.
"""
from __future__ import annotations

__all__ = ["_none_leq", "_none_geq", "_none_eq"]


def _none_leq(a: int | float | None, b: int | float | None) -> bool:
    """``a <= b``, treating a ``None`` operand as unbounded (``True``)."""
    if a is None or b is None:
        return True
    return a <= b


def _none_geq(a: int | float | None, b: int | float | None) -> bool:
    """``a >= b``, treating a ``None`` operand as unbounded (``True``)."""
    if a is None or b is None:
        return True
    return a >= b


def _none_eq(a: int | float | None, b: int | float | None) -> bool:
    """``a == b``, treating a ``None`` operand as a wildcard (``True``)."""
    if a is None or b is None:
        return True
    return a == b
