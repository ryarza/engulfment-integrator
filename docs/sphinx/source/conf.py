# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'engulfment-integrator'
copyright = '2022, Ricardo Yarza'
author = 'Ricardo Yarza'

# The full version, including alpha/beta/rc tags
#release = '1.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.imgmath', 
    'sphinx.ext.todo',
    'sphinx_rtd_theme',
    'breathe',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

html_theme = "sphinx_rtd_theme"
#html_theme = 'alabaster'
html_theme_options = {
    "collapse_navigation" : False,
}
html_logo = "insp.png"

html_context = {
  'display_github': True,
  'github_user': 'ryarza',
  'github_repo': 'engulfment-integrator'
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']
html_static_path = []


import subprocess, os
subprocess.call('make clean', shell=True)
subprocess.call('cd ../../doxygen ; doxygen', shell=True)

breathe_projects = { project : "../../doxygen/build/xml/"}#, "star": "../../doxygen/build/xml/" }

#Find all header files in the code
header_files = []
for file in os.listdir("../../../src"):
    if file.endswith(".h"):
        header_files.append(file)

breathe_default_project = project
breathe_projects_source = {
     project : ( "../../../src/", header_files )
     }
