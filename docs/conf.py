# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'FinaleToolkit'
copyright = '2024, EpiFluidLab'
author = 'EpiFluidLab'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

extensions = [
    'sphinx.ext.autodoc',
    'sphinxarg.ext',
    'sphinx.ext.napoleon'
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_logo = "_static/finaletoolkit_logo_rounded.png"
html_favicon = "_static/favicon.ico"
html_theme_options = {
    "logo": {
        "text": "FinaleToolkit",
        "image_dark": "_static/finaletoolkit_logo_rounded.png",
    }
}

