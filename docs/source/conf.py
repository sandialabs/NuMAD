# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'NuMAD'
copyright = 'Copyright 2013-2021, National Technology & Engineering Solutions of Sandia, LLC (NTESS).'
author = 'NumAD Authors'

release = '3.0'
version = '3.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
