"""
Shared multiprocessing and progress-reporting helpers.

This module centralizes the ``multiprocessing.Pool`` + ``tqdm`` pattern that was
copy-pasted across the original feature modules so that worker-count handling
and progress display are consistent everywhere.
"""
from __future__ import annotations

import multiprocessing as mp
from contextlib import contextmanager
from typing import Callable, Iterable, Iterator, Sequence, TypeVar

from tqdm import tqdm

__all__ = ["resolve_workers", "chunksize_for", "progress", "pool", "imap_unordered_progress"]

_T = TypeVar("_T")
_R = TypeVar("_R")


def resolve_workers(workers: int | None) -> int:
    """Return a sane worker count (>= 1).

    Parameters
    ----------
    workers : int or None
        Requested worker count.  ``None`` or values < 1 fall back to 1.
    """
    if workers is None or workers < 1:
        return 1
    return int(workers)


def chunksize_for(n_items: int, workers: int, cap: int | None = None) -> int:
    """Compute a balanced ``Pool.imap`` chunksize.

    Mirrors the original ``max(len(items) // workers, 1)`` heuristic, with an
    optional upper ``cap``.

    Parameters
    ----------
    n_items : int
        Number of work items.
    workers : int
        Number of worker processes.
    cap : int, optional
        Maximum chunksize.
    """
    size = max(n_items // max(workers, 1), 1)
    if cap is not None:
        size = min(size, cap)
    return size


def progress(
    iterable: Iterable[_T],
    *,
    verbose: bool | int = False,
    desc: str | None = None,
    total: int | None = None,
    **kwargs,
) -> Iterable[_T]:
    """Wrap ``iterable`` in a tqdm progress bar when ``verbose`` is truthy.

    When ``verbose`` is falsy the iterable is returned unchanged, so callers can
    always iterate the result.

    Parameters
    ----------
    iterable : iterable
        The iterable to consume.
    verbose : bool or int, optional
        Show a progress bar when truthy.
    desc : str, optional
        Progress-bar description.
    total : int, optional
        Expected item count (for iterables without ``len``).
    **kwargs
        Forwarded to :class:`tqdm.tqdm`.
    """
    if not verbose:
        return iterable
    return tqdm(iterable, desc=desc, total=total, **kwargs)


@contextmanager
def pool(workers: int, maxtasksperchild: int | None = None) -> Iterator[mp.pool.Pool]:
    """Context manager yielding a closed-and-joined :class:`multiprocessing.Pool`.

    Parameters
    ----------
    workers : int
        Number of worker processes (coerced to >= 1).
    maxtasksperchild : int, optional
        Recycle workers after this many tasks (bounds memory growth).
    """
    p = mp.Pool(processes=resolve_workers(workers), maxtasksperchild=maxtasksperchild)
    try:
        yield p
    finally:
        p.close()
        p.join()


def imap_unordered_progress(
    func: Callable[[_T], _R],
    items: Sequence[_T],
    workers: int,
    *,
    chunksize: int | None = None,
    maxtasksperchild: int | None = None,
    verbose: bool | int = False,
    desc: str | None = None,
) -> list[_R]:
    """Map ``func`` over ``items`` with a pool, returning a list of results.

    Convenience wrapper for the common "fan out, collect, maybe show progress"
    pattern.  Use the lower-level :func:`pool` when result ordering or streaming
    matters.

    Parameters
    ----------
    func : callable
        Worker function applied to each item.
    items : sequence
        Work items.
    workers : int
        Worker-process count.
    chunksize : int, optional
        ``imap`` chunksize; computed from ``items``/``workers`` if omitted.
    maxtasksperchild : int, optional
        Recycle workers after this many tasks.
    verbose : bool or int, optional
        Show a progress bar when truthy.
    desc : str, optional
        Progress-bar description.

    Returns
    -------
    list
        Results in completion order.
    """
    n = len(items)
    cs = chunksize if chunksize is not None else chunksize_for(n, workers)
    with pool(workers, maxtasksperchild=maxtasksperchild) as p:
        return list(
            progress(
                p.imap_unordered(func, items, chunksize=cs),
                verbose=verbose,
                total=n,
                desc=desc,
            )
        )
