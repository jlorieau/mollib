Mollib
======

.. image:: https://travis-ci.org/jlorieau/mollib.svg?branch=master
    :target: https://travis-ci.org/jlorieau/mollib

Mollib is a *unified* command-line program and Python library for the
validation, quality analysis and manipulation of molecular structures with an
emphasis on biophysical analysis. Mollib is built on a plugin framework to
easily add new tools to manipulate and analyze structures and data, which can
then be *combined* and *cross-validated*.

It includes tools for:

- The `processing`_ and protonation of molecules.
- The analysis of `hydrogen bonds`_ and hydrogen bond quality compared to the
  highest-resolution PDB structures. Hydrogen bonds are classified based on
  their acceptor-donor residues and Ramachandran angles
- The `measurement`_ of geometries within molecules as well as the reporting
  and classification of Ramachandran angles.
- The analysis `partial alignment`_ NMR data with residual dipolar coupling
  (RDC) and residual anisotropic chemical shift (RACS, a.k.a RCSA) data.

.. _`processing`: http://mollib.readthedocs.io/en/latest/cli/process.html
.. _`hydrogen bonds`: http://mollib.readthedocs.io/en/latest/cli/hbonds.html
.. _`measurement`: http://mollib.readthedocs.io/en/latest/cli/measure.html
.. _`partial alignment`: http://mollib.readthedocs.io/en/latest/cli/pa.html

Getting Started
---------------

Mollib can be installed through pip.

.. code:: shell

    $ pip install mollib

The mollib documentation (`html`_ | `pdf`_) is generously hosted by readthedocs.


The `mollib source code`_ is hosted on github.

.. _`html`: http://mollib.readthedocs.io/en/latest/
.. _`pdf`: http://readthedocs.org/projects/mollib/downloads/pdf/latest/
.. _`mollib source code`: https://github.com/jlorieau/mollib