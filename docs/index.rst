mollib's documentation
######################


Mollib is a *unified* command-line program and Python library for the
validation, quality analysis and manipulation of molecular structures with an
emphasis on biophysical analysis. Mollib is built on a plugin framework to
easily add new tools to manipulate and analyze structures and data, which can
then be *combined* and *cross-validated*.

It includes tools for:

- The :ref:`processing <process-command>` and protonation of
  molecules.
- The analysis of :ref:`hydrogen bonds <hbonds-command>` and hydrogen bond
  quality compared to the highest-resolution PDB structures. Hydrogen bonds
  are classified based on their acceptor-donor residues and Ramachandran
  angles
- The :ref:`measurement <measure-command>` of geometries within molecules
  as well as the reporting and classification of Ramachandran angles.
- The analysis :ref:`partial alignment <pa-command>` data with residual dipolar
  coupling (RDC) and residual anisotropic chemical shift (RACS, a.k.a RCSA)
  data.

Installation
============

**Github**: https://github.com/jlorieau/mollib

**Download**: `mollib-1.3.zip <https://github.com/jlorieau/mollib/archive/v1.3.zip>`_
`mollib-1.3.tar.gz <https://github.com/jlorieau/mollib/archive/v1.3.tar.gz>`_



Indices and tables
==================

.. toctree::
   :maxdepth: 2

   cli/cli
   api/api
   develop
   releases/releases

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

