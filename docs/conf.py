# Configuration file for the Sphinx documentation builder.
#
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Make the package importable for autodoc even when it is not installed
# (src layout). When installed in the build environment this is redundant but
# harmless.
sys.path.insert(0, os.path.abspath("../src"))

# -- Project information -----------------------------------------------------

project = "FinaleToolkit"
copyright = "2026, EpiFluidLab"
author = "EpiFluidLab"

# -- General configuration ---------------------------------------------------

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinxarg.ext",
]

# NumPy-style docstrings.
napoleon_numpy_docstring = True
napoleon_google_docstring = False

# Autodoc: keep source order, show type hints in the description.
autodoc_member_order = "bysource"
autodoc_typehints = "description"

# Note: autodoc imports the package, so the build environment must have the
# runtime dependencies (numpy, pysam, numba, ...) installed. They are not
# mocked because numba's @jit decorators need the real library to import.

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
}

# -- Options for HTML output -------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_logo = "_static/finaletoolkit_logo_rounded.png"
html_favicon = "_static/favicon.ico"
html_theme_options = {
    "logo": {
        "text": "FinaleToolkit",
        "image_dark": "_static/finaletoolkit_logo_rounded.png",
    }
}
