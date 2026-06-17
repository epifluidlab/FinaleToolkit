"""
Lightweight logging wrapper used across FinaleToolkit.

A single stderr handler is attached per named logger with a consistent
``[timestamp] LEVEL [name] message`` format.  ``set_verbosity`` adjusts the
shared ``finaletoolkit`` parent logger so child loggers inherit the level.
"""
from __future__ import annotations

import logging
import sys
from typing import Any

__all__ = ["Logger", "get_logger", "set_verbosity"]


class Logger:
    """A thin, consistent wrapper around :class:`logging.Logger`."""

    def __init__(self, name: str, level: int = logging.INFO) -> None:
        self._logger = logging.getLogger(name)
        self._setup_handler(level)

    def _setup_handler(self, level: int) -> None:
        """Attach a stderr handler with package formatting (once per logger)."""
        if not self._logger.handlers:
            self._logger.setLevel(level)

            handler = logging.StreamHandler(sys.stderr)
            handler.setLevel(level)
            formatter = logging.Formatter(
                fmt="[%(asctime)s] %(levelname)s [%(name)s] %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
            handler.setFormatter(formatter)

            self._logger.addHandler(handler)
            self._logger.propagate = False

    def debug(self, msg: Any, *args: Any, **kwargs: Any) -> None:
        """Log a message with severity ``DEBUG``."""
        self._logger.debug(msg, *args, **kwargs)

    def info(self, msg: Any, *args: Any, **kwargs: Any) -> None:
        """Log a message with severity ``INFO``."""
        self._logger.info(msg, *args, **kwargs)

    def warning(self, msg: Any, *args: Any, **kwargs: Any) -> None:
        """Log a message with severity ``WARNING``."""
        self._logger.warning(msg, *args, **kwargs)

    def error(self, msg: Any, *args: Any, **kwargs: Any) -> None:
        """Log a message with severity ``ERROR``."""
        self._logger.error(msg, *args, **kwargs)

    def critical(self, msg: Any, *args: Any, **kwargs: Any) -> None:
        """Log a message with severity ``CRITICAL``."""
        self._logger.critical(msg, *args, **kwargs)

    def set_level(self, level: int) -> None:
        """Set the logging level for this logger and its handlers."""
        self._logger.setLevel(level)
        for handler in self._logger.handlers:
            handler.setLevel(level)


def get_logger(name: str, level: int = logging.INFO) -> Logger:
    """Return a :class:`Logger` for ``name`` (typically ``__name__``).

    Parameters
    ----------
    name : str
        The logger name.
    level : int, optional
        The logging level, by default :data:`logging.INFO`.

    Returns
    -------
    Logger
        The wrapped logger instance.
    """
    return Logger(name, level)


def set_verbosity(level: int) -> None:
    """Set the level of the shared ``finaletoolkit`` parent logger.

    All child loggers inherit this level unless explicitly overridden.

    Parameters
    ----------
    level : int
        A logging level such as :data:`logging.DEBUG` or :data:`logging.INFO`.
    """
    logging.getLogger("finaletoolkit").setLevel(level)
