######################
mollib's documentation
######################


MolLib is a Python module and commandline program for the validation, quality
analysis and manipulation molecular structures with an emphasis on biophysical
analysis. The main objective are:

- A *powerful* and *unified* Python package for analyzing and
  manipulating structures.
- A *powerful* and *unified* commandline program for analyzing and manipulating
  structures.
- To make a plugin framework to easily add tools for the validation, analysis
  and manipulation of structures.
- By using a common framework, tools can be *cross-validated* and *combined*.
- All in one place: Provide standard and new quality analysis and
  manipulation methods provided by Molprobity, Procheck and other
  tools into one package.
- Provide an open source reference implementation of analysis and manipulation
  methods to document and teach how these methods are conducted.
- Provide common functions and object structures for molecular analyses that
  are computationally efficient (fast) and parallelizable. Computational speed
  is achieved by converting time consuming portions of the program into C, C++
  and Cython.
- Produce detailed, easy to read reports on structures

**Github**: https://github.com/jlorieau/mollib

**Download**: `mollib-1.3.zip <https://github.com/jlorieau/mollib/archive/v1.3.zip>`_
`mollib-1.3.tar.gz <https://github.com/jlorieau/mollib/archive/v1.3.tar.gz>`_

Features
========

* Plugin interface to easily add functionality.
* Common manipulations of biomolecules for biophysics, like translation,
  Euler rotations, hydrogenation and hydrogen bond detection.
* Easily extendable data types (Atom, Residue, Chain, Molecule) through object
  inheritance to add functionality and data to base objects.
* Loading, fetching, caching  and writing PDB files.
* Handles gzipped and and plain text files
* A powerful commandline interface and python module library

Description
===========

.. automodule:: mollib.core.molecule

Indices and tables
==================

.. toctree::
   :maxdepth: 3

   cli/cli
   api/api
   develop


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

