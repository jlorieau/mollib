
# -*- coding: utf-8 -*-
#
# MolLib documentation build configuration file, created by
# sphinx-quickstart on Thu Aug  4 18:54:31 2016.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import textwrap
sys.path.insert(0, os.path.abspath('..'))

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinxcontrib.napoleon',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The encoding of source files.
#
# source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = 'mollib'
copyright = u'2016, Justin L Lorieau'
author = 'Justin L Lorieau'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.

# Get the version number and package information.
# The __version__.py file is executed so that the mollib package is not loaded.
# At this point, the C/C++ extensions may not be built, and loading the mollib
# package will lead to an ImportError. This approach circumvents this problem.
__version__ = None  # This is a version string
VERSION = None  # This is a 5-item version tuple
exec(open("../mollib/__version__.py").read())

# The full version, including alpha/beta/rc tags.
release = __version__


# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
# language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#
# today = ''
#
# Else, today_fmt is used as the format for a strftime call.
#
# today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'cli/cmds']

# The reST default role (used for this markup: `text`) to use for all
# documents.
#
# default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#
# add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#
# add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#
# show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# If true, keep warnings as "system message" paragraphs in the built documents.
# keep_warnings = False

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
#html_theme = 'sphinx_rtd_theme'
html_theme = 'nature'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {'sidebarwidth': 200}

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = []

# The name for this set of Sphinx documents.
# "<project> v<release> documentation" by default.
#
html_title = u'Mollib ' + __version__

# A shorter title for the navigation bar.  Default is the same as html_title.
#
# html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#
# html_logo = None

# The name of an image file (relative to this directory) to use as a favicon of
# the docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#
# html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Add any extra paths that contain custom files (such as robots.txt or
# .htaccess) here, relative to this directory. These files are copied
# directly to the root of the documentation.
#
# html_extra_path = []

# If not None, a 'Last updated on:' timestamp is inserted at every page
# bottom, using the given strftime format.
# The empty string is equivalent to '%b %d, %Y'.
#
# html_last_updated_fmt = None

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#
# html_sidebars = {}
# html_sidebars = { '**': ['globaltoc.html', 'relations.html', 'sourcelink.html',
#                          'searchbox.html'], }
html_sidebars = { '**': ['globaltoc.html', ], }

# Additional templates that should be rendered to pages, maps page names to
# template names.
#
# html_additional_pages = {}

# If false, no module index is generated.
#
html_domain_indices = True

# If false, no index is generated.
#
html_use_index = True

# If true, the index is split into individual pages for each letter.
#
# html_split_index = False

# If true, links to the reST sources are added to the pages.
#
# html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#
# html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#
html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#
# html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = None

# Language to be used for generating the HTML full-text search index.
# Sphinx supports the following languages:
#   'da', 'de', 'en', 'es', 'fi', 'fr', 'hu', 'it', 'ja'
#   'nl', 'no', 'pt', 'ro', 'ru', 'sv', 'tr', 'zh'
#
# html_search_language = 'en'

# A dictionary with options for the search language support, empty by default.
# 'ja' uses this config value.
# 'zh' user can custom change `jieba` dictionary path.
#
# html_search_options = {'type': 'default'}

# The name of a javascript file (relative to the configuration directory) that
# implements a search results scorer. If empty, the default will be used.
#
# html_search_scorer = 'scorer.js'

# Output file base name for HTML help builder.
htmlhelp_basename = 'mollibdoc'

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
     # The paper size ('letterpaper' or 'a4paper').
     #
     # 'papersize': 'letterpaper',

     'classoptions': ',openany,oneside',

     'fontpkg': r'''
        \usepackage{times}
        \usepackage{inconsolata}
      ''',

     # The font size ('10pt', '11pt' or '12pt').
     #
     # 'pointsize': '11pt',

     # Additional stuff for the LaTeX preamble.
     #
     'preamble': r'''
         \definecolor{olivegreen}{RGB}{60, 128, 49}
         \definecolor{darkyellow}{RGB}{202, 214, 41}
         \definecolor{darkorange}{RGB}{198,  93,  9}
         \sphinxDeclareColorOption{VerbatimBorderColor}{{rgb}{0.5,0.5,0.5}}
         \sphinxDeclareColorOption{VerbatimColor}{{rgb}{0.985,0.985,0.985}}
         \fvset{fontsize=auto}
         \sphinxverbatimsep=6pt
         \sphinxshadowsize=15pt
         \usepackage{enumitem}
         \setlistdepth{99}
         \setlist{labelsep=0.5em}
         \DeclareUnicodeCharacter{2212}{-}
         '''

     # Latex figure (float) alignment
     #
     # 'figure_align': 'htbp',

}

# The API documentation is not included in the user PDF manual
if 'tags' in locals() and tags.has('latex'):
    exclude_patterns += ['**/releases*',]
    # exclude_patterns += ['**/api*', '**/releases*', 'develop*']

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    ('index', 'mollib.tex',
     u'Mollib Documentation',
     u'Justin L Lorieau', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#
# latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#
latex_use_parts = False

# If true, show page references after internal links.
#
# latex_show_pagerefs = False

# If true, show URL addresses after external links.
#
# latex_show_urls = False

# Documents to append as an appendix to all manuals.
#
# latex_appendices = []

# It false, will not define \strong, \code, 	itleref, \crossref ... but only
# \sphinxstrong, ..., \sphinxtitleref, ... To help avoid clash with user added
# packages.
#
# latex_keep_old_macro_names = True

# If false, no module index is generated.
#
# latex_domain_indices = True


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'mollib', u'Mollib Documentation',
     [author], 1)
]

# If true, show URL addresses after external links.
#
# man_show_urls = False


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'MolLib', u'MolLib Documentation',
     author, 'MolLib', 'One line description of project.',
     'Miscellaneous'),
]

# Documents to append as an appendix to all manuals.
#
# texinfo_appendices = []

# If false, no module index is generated.
#
# texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#
# texinfo_show_urls = 'footnote'

# If true, do not generate a @detailmenu in the "Top" node's menu.
#
# texinfo_no_detailmenu = False

# Sort autodoc elements by source, not alphabetically (default)
autodoc_member_order = 'bysource'

# Custom css
def setup(app):
    app.add_stylesheet('updates.css')

# Get the CLI output text
def process_cmd(string):
    # Process the command. Split the arguments and strip '-' and '|' characters
    args = string.split()
    progname = args[0]
    args_name = '_'.join([i.strip('-').replace('|', '_') for i in args[1:]])
    shell_cmd = "$ {cmd}".format(cmd=string)
    shell_cmd = " \\\n> ".join(textwrap.wrap(shell_cmd, 78))

    # Prepare the CLI output file.
    print(shell_cmd)
    cmd = ("echo '{shell_cmd}' "
            "> cli/output/{progname}_{args_name}.txt\n")

    cmd += ("cd ..&&"
            "FORCE_COLOR=TRUE {progname} {args}"
            ">> docs/cli/output/{progname}_{args_name}.txt\n"
            "cd docs\n")

    # Process the html component
    cmd += ("echo '.. only:: html\n\n.. raw:: html\n'"
            "> cli/output/{progname}_{args_name}.html\n")

    cmd += ("pygmentize -l shell-session -f html "
            "cli/output/{progname}_{args_name}.txt"
            "|sed 's/^/    /g' >> cli/output/{progname}_{args_name}.html\n")

    # Replace ANSI colors
    cmd += ("cat -e cli/output/{progname}_{args_name}.html"
            "|sed 's/\$$//g'"  # Remove $ at the end of lines
            "|sed 's/\^\[\[1m/<font style=\"font-weight:bold;\">/g'"
            "|sed 's/\^\[\[22m/<\/font>/g'"
            "|sed 's/\^\[\[36m/<font color=\"#008b8b\">/g'"  # cyan
            "|sed 's/\^\[\[91m/<font color=\"red\">/g'"
            "|sed 's/\^\[\[92m/<font color=\"green\">/g'"
            "|sed 's/\^\[\[33m/<font color=\"#abb51f\">/g'"
            "|sed 's/\^\[\[94m/<font color=\"blue\">/g'"
            "|sed 's/\^\[\[95m/<font color=\"magenta\">/g'"
            "|sed 's/\^\[\[96m/<font color=\"#008b8b\">/g'"  # cyan
            "|sed 's/\^\[\[0m/<\/font>/g'"
            ">cli/output/{progname}_{args_name}.tmp\n")

    cmd += ("mv cli/output/{progname}_{args_name}.tmp "
            "cli/output/{progname}_{args_name}.rst\n")

    # Process the latex component
    cmd += ("echo '\n.. only:: latex\n\n"
            ".. raw:: latex\n\n"
            "    \\\\begin{{sphinxVerbatim}}"
            "[commandchars=\\\\\\\\\\\\{{\\\\}},"  # [commandchars=\\\{\},
            "fontsize=\\\\footnotesize]"  # fontsize=\footnotesize,
            "'"
            "> cli/output/{progname}_{args_name}.tex\n")

    cmd += ("cat -e cli/output/{progname}_{args_name}.txt"
            "|sed 's/user@host\$/\\\\textcolor{{darkorange}}"
                "{{\\\\textbf{{user@host$}}}}/g'"  # highlight the term prompt
            "|sed 's/^/    /g'"  # Add a tab at the start of every line
            "|sed 's/^   //g'" # Strip leading spaces
            "|sed 's/\$$//g'"  # Remove $ at the end of lines
            "|sed 's/--/-{{-}}/g'"  # Preserve --
            "|sed 's/\^\[\[1m/\\\\textbf{{/g'"
            "|sed 's/\^\[\[22m/}}/g'"
            "|sed 's/\^\[\[36m/\\\\textcolor{{cyan}}{{/g'"
            "|sed 's/\^\[\[91m/\\\\textcolor{{red}}{{/g'"
            "|sed 's/\^\[\[92m/\\\\textcolor{{olivegreen}}{{/g'"
            "|sed 's/\^\[\[33m/\\\\textcolor{{darkyellow}}{{/g'"
            "|sed 's/\^\[\[94m/\\\\textcolor{{blue}}{{/g'"
            "|sed 's/\^\[\[95m/\\\\textcolor{{magenta}}{{/g'"
            "|sed 's/\^\[\[96m/\\\\textcolor{{cyan}}{{/g'"
            "|sed 's/\^\[\[0m/}}/g'"
            ">>cli/output/{progname}_{args_name}.tex\n")

    #cmd += ("echo '\n    \\\\end{{sphinxVerbatim}}\n\\\\\\\\\n'"
    cmd += ("echo '    \\\\end{{sphinxVerbatim}}\n {{}} \n'"
            ">>cli/output/{progname}_{args_name}.tex\n")

    cmd += ("cat cli/output/{progname}_{args_name}.tex "
            ">> cli/output/{progname}_{args_name}.rst\n")

    # Clean up
    cmd += ("rm cli/output/*.txt cli/output/*.tex cli/output/*.html\n")
    cmd = cmd.format(shell_cmd=shell_cmd, progname=progname,
                     args=' '.join(args[1:]), args_name=args_name)

    os.system(cmd)

if 'cli' in sys.argv:
    process_cmd("make help")
    #
    # process_cmd("ml --help")
    # process_cmd("ml --list-plugins")
    # process_cmd("ml --list-settings")
    #
    # process_cmd("ml process --help")
    # process_cmd("ml process -i 1UBQ -o 1UBQ_H.pdb --hydrogenate")

    # process_cmd("ml measure --help")
    # process_cmd("ml measure -i 2MUV -d A:D.20:21.CA A:D.20:21.CA --exclude-intra --exclude-intra-chain")
    # process_cmd("ml measure -i 2MUV -d 23:49.HA 23:49.H --only-delta 3 --stats")
    # process_cmd("ml measure -i 2KXA 2LWA -d A:C.5.HA A:C.21.H --only-intra-chain")
    # process_cmd("ml measure -i 2MJB -a 20:30.C 20:30.N 20:30.H --only-bonded")
    # process_cmd("ml measure -i 2KXA -dih 2:6.C 2:6.N 2:6.CA 2:6.C --only-bonded --stats")
    # process_cmd("ml measure -i 2KXA --rama")

    # process_cmd("ml hbonds --help")
    # process_cmd("ml hbonds -i 2KXA")
    # process_cmd("ml hbonds -i 1UBQ --hydrogenate|head -n15")

    # process_cmd("ml pa --help")
    # process_cmd("ml pa -i 2KXA -a 2KXA")
    # process_cmd("ml pa -i 2MJB -a 2MJB --set 0 --fix-outliers --project-methyls --summary")
    # process_cmd("ml pa -i 2MJB -a 2MJB --set 0 --exclude CE-HE CD-HD CE-SD --fix-outliers --project-methyls --summary")
    # process_cmd("ml pa -i 1UBQ -a 2MJB --set 0 --fix-outliers --project-methyls --hydrogenate --summary")

    # cleanup
    os.system("rm ../1UBQ_H.pdb")
