# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'KGGSEE User Manual'
copyright = '2019 PMGLab'
author = 'Miaoxin Li, Lin Jiang, Xiangyi Li, and Lin Miao'

# release = ''
version = 'v1.1'

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

# -- Options for LaTeX output
latex_elements = {
    'extraclassoptions': 'openany,oneside',
    'papersize': 'a4paper',
    'pointsize': '10pt',
}

html_static_path = ['_static']
def setup(app):
    app.add_css_file('css/kggsee.css')
