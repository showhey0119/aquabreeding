# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

project = 'aquabreeding'
copyright = '2022, Shohei Takuno'
author = 'Shohei Takuno'
release = '0.8.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx_rtd_theme'
    ]

[extensions]
todo_include_todos=True

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
# html_theme = 'sphinx_rtd_theme'
# html_theme = 'nature'
html_static_path = ['_static']

html_theme_options = {
    "light_css_variables": {
        "color-admonition-background": "#E5FFE5",  # background of Note, Todo
        "color-foreground-primary": "#001100",  # main text
        "color-foreground-muted": "#006600",  # contents
        "color-background-hover": "#CCFFCC",  # hover
        "color-background-border": "#007700",  # object, underline of link
        "color-background-item": "#FF0000",  # underline of link
        "color-problematic": "#007700",  # class color
        "api-font-size": "var(--font-size--normal)",  # class font size
    },
}

# use this to modify css
# html_style = "css/my_theme.css"


