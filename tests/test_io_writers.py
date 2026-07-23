"""
Tests for finaletoolkit.io.writers
"""

import gzip
import sys

from finaletoolkit.io.writers import smart_open_text, is_stdout


class TestIsStdout:
    def test_dash_is_stdout(self):
        assert is_stdout("-")

    def test_path_is_not_stdout(self):
        assert not is_stdout("output.txt")


class TestSmartOpenText:
    def test_writes_stdout(self, capsys):
        with smart_open_text("-") as f:
            assert f is sys.stdout
            f.write("hello\n")
        assert capsys.readouterr().out == "hello\n"

    def test_stdout_not_closed_on_exit(self, capsys):
        with smart_open_text("-") as f:
            pass
        assert not sys.stdout.closed

    def test_writes_plain_text_file(self, tmp_path):
        path = tmp_path / "out.txt"
        with smart_open_text(str(path)) as f:
            f.write("plain text\n")
        assert path.read_text() == "plain text\n"

    def test_writes_gzip_file(self, tmp_path):
        path = tmp_path / "out.txt.gz"
        with smart_open_text(str(path)) as f:
            f.write("gzipped text\n")
        with gzip.open(path, "rt") as f:
            assert f.read() == "gzipped text\n"

    def test_file_closed_on_exit(self, tmp_path):
        path = tmp_path / "out.txt"
        with smart_open_text(str(path)) as f:
            handle = f
        assert handle.closed

    def test_file_closed_on_exception(self, tmp_path):
        path = tmp_path / "out.txt"
        handle = None
        try:
            with smart_open_text(str(path)) as f:
                handle = f
                raise ValueError("boom")
        except ValueError:
            pass
        assert handle.closed
