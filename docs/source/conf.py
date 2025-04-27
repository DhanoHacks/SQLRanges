import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SQLRanges'
copyright = '2025, Dhananjay Raman, Saket Choudhary'
author = 'Dhananjay Raman, Saket Choudhary'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',  # for Google-style or NumPy-style docstrings
    'sphinx.ext.viewcode',  # (optional) adds "view source" links
]

templates_path = ['_templates']
exclude_patterns = []

autodoc_default_flags = ['members', 'undoc-members', 'show-inheritance']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_baseurl = "https://dhanohacks.github.io/SQLRanges/"