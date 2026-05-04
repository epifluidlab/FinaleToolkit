import logging
import sys
from typing import Optional, Any

class Logger:
    """
    A professional wrapper around the standard Python logging.Logger.
    Provides a cleaner interface and consistent formatting across the toolkit.
    """
    def __init__(self, name: str, level: int = logging.INFO):
        self._logger = logging.getLogger(name)
        self._setup_handler(level)

    def _setup_handler(self, level: int):
        """Sets up the professional formatting if not already configured."""
        if not self._logger.handlers:
            self._logger.setLevel(level)
            
            # Create console handler
            handler = logging.StreamHandler(sys.stderr)
            handler.setLevel(level)
            
            # Professional format: [TIMESTAMP] LEVEL [MODULE] MESSAGE
            formatter = logging.Formatter(
                fmt='[%(asctime)s] %(levelname)s [%(name)s] %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
            handler.setFormatter(formatter)
            
            self._logger.addHandler(handler)
            self._logger.propagate = False

    def debug(self, msg: Any, *args: Any, **kwargs: Any):
        """Log a message with severity 'DEBUG'."""
        self._logger.debug(msg, *args, **kwargs)

    def info(self, msg: Any, *args: Any, **kwargs: Any):
        """Log a message with severity 'INFO'."""
        self._logger.info(msg, *args, **kwargs)

    def warning(self, msg: Any, *args: Any, **kwargs: Any):
        """Log a message with severity 'WARNING'."""
        self._logger.warning(msg, *args, **kwargs)

    def error(self, msg: Any, *args: Any, **kwargs: Any):
        """Log a message with severity 'ERROR'."""
        self._logger.error(msg, *args, **kwargs)

    def critical(self, msg: Any, *args: Any, **kwargs: Any):
        """Log a message with severity 'CRITICAL'."""
        self._logger.critical(msg, *args, **kwargs)

    def set_level(self, level: int):
        """Set the logging level for this specific logger."""
        self._logger.setLevel(level)
        for handler in self._logger.handlers:
            handler.setLevel(level)

def get_logger(name: str, level: int = logging.INFO) -> Logger:
    """
    Returns a Logger instance with the specified name and level.
    
    Parameters
    ----------
    name : str
        The name of the logger, typically __name__.
    level : int, optional
        The logging level, by default logging.INFO.
        
    Returns
    -------
    Logger
        The wrapped logger instance.
    """
    return Logger(name, level)

def set_verbosity(level: int):
    """
    Sets the verbosity level for the 'finaletoolkit' parent logger.
    All child loggers will inherit this level unless specifically set.
    
    Parameters
    ----------
    level : int
        The logging level (e.g., logging.DEBUG, logging.INFO).
    """
    logging.getLogger('finaletoolkit').setLevel(level)
