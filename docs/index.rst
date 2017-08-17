######
mollib
######


Mollib is a *unified* command-line program and Python library for the
validation, quality analysis and manipulation of molecular structures with an
emphasis on biophysical analysis. Mollib is built on a plugin framework to
easily add new tools to manipulate and analyze structures and data, which can
then be *combined* and *cross-validated*.

Mollib includes tools for:

- The :ref:`processing <process-command>` and protonation of
  molecules.
- The analysis of :ref:`hydrogen bonds <hbonds-command>` and hydrogen bond
  quality compared to the highest-resolution PDB structures. Hydrogen bonds
  are classified based on their acceptor-donor residues and Ramachandran
  angles
- The :ref:`measurement <measure-command>` of geometries within molecules
  as well as the reporting and classification of Ramachandran angles.
- The statistical analysis and comparison of structures to high-resolution
  crystal structures.
- The analysis of :ref:`partial alignment <pa-command>` data with residual
  dipolar coupling (RDC) and residual anisotropic chemical shift (RACS, a.k.a
  RCSA) data.


.. rubric:: Installation

.. code:: shell

    $ pip install mollib

.. note:: Python and pip may not be installed. You may need to
          `download <https://www.python.org/downloads/>`_ and install
          Python first.

.. only:: html

    .. rubric:: Table of Contents

.. toctree::
   :maxdepth: 2

   cli/cli
   recipes/recipes
   releases/releases
   develop
   building
   api/api
