import warnings
from functools import wraps


def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used. A similar decorator was introduced in
    Python 3.13, which if backported should be used instead."""
    @wraps(func)
    def new_func(*args, **kwargs):
        warnings.warn(f"Call to deprecated function {func.__name__}.",
                      category=DeprecationWarning,
                      stacklevel=2)
        return func(*args, **kwargs)
    return new_func


def moved(new_function):
    """
    A decorator to mark a function as renamed and deprecated.
    It warns the user and redirects calls to the new function.

    Parameters
    ==========
    new_function: callable
        The new function that replaces the deprecated one.
    """
    def decorator(old_function):
        @wraps(old_function)
        def wrapped(*args, **kwargs):
            warnings.warn(
                f"{old_function.__name__} is deprecated and has been renamed "
                f"to {new_function.__name__}. Please use the new function "
                "instead.",
                DeprecationWarning,
                stacklevel=2
            )
            return new_function(*args, **kwargs)
        return wrapped
    return decorator
