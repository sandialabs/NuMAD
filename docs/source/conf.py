# File: docs/source/conf.py 
# Configuration file for the Sphinx documentation builder. 

import os
import re

# -- Project information

project = 'NuMAD'
copyright = '2013-2021, National Technology & Engineering Solutions of Sandia, LLC (NTESS)' 
author = 'The NuMAD Developers'

version = 'v3.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinxcontrib.bibtex',
    'sphinxcontrib.matlab',
    'sphinxext.remoteliteralinclude',       
]

#intersphinx_mapping = {
#    'python': ('https://docs.python.org/3/', None),
#    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
#}
#intersphinx_disabled_domains = ['std']

# autodoc settings
autodoc_member_order = 'bysource'

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# sphinxcontrib.bibtex settings
bibtex_bibfiles = ['refs/userGuide.bib']

# sphinxcontrib.matlab settings
primary_domain = 'mat'
matlab_src_dir = os.path.abspath("../../source")
matlab_keep_package_prefix = False

# Add any paths that contain templates here, relative to this directory.
#templates_path = ['_templates']

# The suffix(es) of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# Enable numref
numfig = True

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_context = {
    'css_files': [
        '_static/theme_overrides.css',  # override wide tables in RTD theme
        ],
     }


# -- Options for EPUB output
#epub_show_urls = 'footnote'
