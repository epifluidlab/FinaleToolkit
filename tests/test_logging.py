"""
Tests for finaletoolkit.utils.logging
"""

import logging

from finaletoolkit.utils.logging import Logger, get_logger, set_verbosity


class TestLogger:
    def test_get_logger_returns_logger(self):
        log = get_logger("finaletoolkit.test_logging.a")
        assert isinstance(log, Logger)

    def test_handler_attached_once(self):
        # Constructing two Logger wrappers for the same name shouldn't
        # duplicate stderr handlers.
        name = "finaletoolkit.test_logging.b"
        first = Logger(name)
        second = Logger(name)
        assert len(first._logger.handlers) == 1
        assert first._logger is second._logger

    def test_log_levels_write_to_stderr(self, capsys):
        # propagate=False keeps records off pytest's caplog handler, so
        # assert on the actual stderr stream the custom handler writes to.
        log = get_logger("finaletoolkit.test_logging.c", level=logging.DEBUG)
        log.debug("debug msg")
        log.info("info msg")
        log.warning("warning msg")
        log.error("error msg")
        log.critical("critical msg")
        err = capsys.readouterr().err
        for msg in ("debug msg", "info msg", "warning msg", "error msg", "critical msg"):
            assert msg in err

    def test_default_level_filters_debug(self, capsys):
        log = get_logger("finaletoolkit.test_logging.c2")
        log.debug("should not appear")
        log.info("should appear")
        err = capsys.readouterr().err
        assert "should not appear" not in err
        assert "should appear" in err

    def test_set_level_updates_logger_and_handlers(self):
        log = get_logger("finaletoolkit.test_logging.d")
        log.set_level(logging.ERROR)
        assert log._logger.level == logging.ERROR
        for handler in log._logger.handlers:
            assert handler.level == logging.ERROR

    def test_set_verbosity_sets_parent_logger_level(self):
        set_verbosity(logging.WARNING)
        assert logging.getLogger("finaletoolkit").level == logging.WARNING
        # child loggers with no explicit level inherit from the parent
        set_verbosity(logging.INFO)
        assert logging.getLogger("finaletoolkit").level == logging.INFO
