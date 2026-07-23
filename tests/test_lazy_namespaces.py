"""
Tests for the lazy attribute-resolution namespaces:
``finaletoolkit`` (flat public API + submodules + aliases) and
``finaletoolkit.cli`` (lazy ``main_cli`` import), both PEP 562 ``__getattr__``
shims that keep ``import finaletoolkit`` cheap.
"""

import pytest

import finaletoolkit


class TestFlatNamespace:
    def test_submodule_lazy_import(self):
        import finaletoolkit.frag as frag_direct

        assert finaletoolkit.frag is frag_direct

    def test_flat_export_resolves_to_real_function(self):
        from finaletoolkit.frag import coverage as coverage_direct

        assert finaletoolkit.coverage is coverage_direct

    def test_flat_export_from_different_submodule(self):
        from finaletoolkit.genome import GenomeGaps as GenomeGaps_direct

        assert finaletoolkit.GenomeGaps is GenomeGaps_direct

    def test_alias_resolves_to_same_object_as_full_name(self):
        assert finaletoolkit.end_motif is finaletoolkit.end_motifs

    def test_breakpoint_alias(self):
        assert finaletoolkit.breakpoint_motif is finaletoolkit.breakpoint_motifs

    def test_unknown_attribute_raises(self):
        with pytest.raises(AttributeError, match="no attribute"):
            finaletoolkit.not_a_real_symbol

    def test_dir_includes_flat_exports_and_submodules_and_aliases(self):
        names = dir(finaletoolkit)
        assert "coverage" in names
        assert "frag" in names
        assert "end_motif" in names
        assert "__version__" in names


class TestCliLazyImport:
    def test_main_cli_lazy_import(self):
        # Calling __getattr__ directly (rather than via plain attribute
        # access, e.g. `cli.main_cli`) is deliberate: `main_cli` is both the
        # submodule's name and the Click group defined inside it, and once
        # anything else in the process imports the `main_cli` submodule
        # (as other test modules do), Python's ordinary import machinery
        # binds the *submodule* onto `finaletoolkit.cli.main_cli`, shadowing
        # whatever `__getattr__` would return -- exactly the ambiguity the
        # module's docstring calls out. Invoking the function isolates the
        # lazy-resolution logic itself from that import-order dependent
        # shadowing.
        import finaletoolkit.cli as cli
        from finaletoolkit.cli.main_cli import main_cli as main_cli_direct

        assert cli.__getattr__("main_cli") is main_cli_direct

    def test_unknown_attribute_raises(self):
        import finaletoolkit.cli as cli

        with pytest.raises(AttributeError, match="no attribute"):
            cli.__getattr__("not_a_real_symbol")
