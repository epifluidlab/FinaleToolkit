from __future__ import annotations

# None compatible comparison operators
def _none_leq(a: int | float | None, b: int | float | None) -> bool:
    """
    Less than or equals that evaluates True if any argument is None
    """
    if a is None or b is None:
        return True
    else:
        return a <= b


def _none_geq(a: int | float | None, b: int | float | None) -> bool:
    """
    Greater than or equals that evaluates True if any argument is None
    """
    if a is None or b is None:
        return True
    else:
        return a >= b


def _none_eq(a: int | float | None, b: int | float | None) -> bool:
    """
    Equals that evaluates True if any argument is None
    """
    if a is None or b is None:
        return True
    else:
        return a == b
