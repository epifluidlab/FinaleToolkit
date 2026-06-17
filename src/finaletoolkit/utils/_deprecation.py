"""
Decorators and helpers for marking deprecated/renamed API.
"""
from __future__ import annotations

import warnings
from functools import wraps
from typing import Callable, TypeVar

__all__ = ["deprecated", "moved", "resolve_length_aliases"]

_T = TypeVar("_T")


def deprecated(func: Callable[..., _T]) -> Callable[..., _T]:
    """Mark a function as deprecated, emitting a :class:`DeprecationWarning`.

    A similar decorator was added to the standard library in Python 3.13
    (:func:`warnings.deprecated`); prefer that once it is the minimum
    supported version.
    """

    @wraps(func)
    def new_func(*args, **kwargs):
        warnings.warn(
            f"Call to deprecated function {func.__name__}.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        return func(*args, **kwargs)

    return new_func


def moved(new_function: Callable[..., _T]) -> Callable[[Callable], Callable[..., _T]]:
    """Mark a function as renamed; warn and redirect calls to ``new_function``.

    Parameters
    ----------
    new_function : callable
        The replacement function.
    """

    def decorator(old_function: Callable) -> Callable[..., _T]:
        @wraps(old_function)
        def wrapped(*args, **kwargs):
            warnings.warn(
                f"{old_function.__name__} is deprecated and has been renamed "
                f"to {new_function.__name__}. Please use the new function "
                "instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            return new_function(*args, **kwargs)

        return wrapped

    return decorator


def resolve_length_aliases(
    min_length: int | None,
    max_length: int | None,
    fraction_low: int | None,
    fraction_high: int | None,
) -> tuple[int | None, int | None]:
    """Resolve the deprecated ``fraction_low``/``fraction_high`` aliases.

    The original toolkit accepted ``fraction_low``/``fraction_high`` as aliases
    for ``min_length``/``max_length``.  This helper centralizes the handling so
    every feature behaves identically: supplying a deprecated alias emits a
    :class:`DeprecationWarning`, and supplying both the old and new spelling of
    the same bound raises :class:`ValueError`.

    Parameters
    ----------
    min_length, max_length : int or None
        The modern length-bound arguments.
    fraction_low, fraction_high : int or None
        The deprecated aliases.

    Returns
    -------
    tuple of (int or None, int or None)
        The resolved ``(min_length, max_length)``.
    """
    if fraction_low is not None:
        if min_length is not None and min_length != fraction_low:
            raise ValueError(
                "fraction_low (deprecated) and min_length were both specified "
                "with different values. Use min_length only."
            )
        warnings.warn(
            "fraction_low is deprecated. Use min_length instead.",
            DeprecationWarning,
            stacklevel=3,
        )
        min_length = fraction_low

    if fraction_high is not None:
        if max_length is not None and max_length != fraction_high:
            raise ValueError(
                "fraction_high (deprecated) and max_length were both specified "
                "with different values. Use max_length only."
            )
        warnings.warn(
            "fraction_high is deprecated. Use max_length instead.",
            DeprecationWarning,
            stacklevel=3,
        )
        max_length = fraction_high

    return min_length, max_length
