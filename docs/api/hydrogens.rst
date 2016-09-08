****************
Hydrogens Module
****************

.. _hydrogenation-implementation:

Implementation
--------------

The hydrogens plugin includes functions to add hydrogens to molecules using
idealized geometry.

The current implementation has the following features:

- adds hydrogens in an idealized geometry with respect to local geometry
  groups.
- Ionizeable groups, like Asp and Glu carboxylates, are protonated based on the
  molecule's sample pH and the set pKa values. When more than one heavy atom
  can be protonated, the heavy atom closest to a hydrogen bond accector will
  be selected.
- Molecules with heteroatoms are correctly hydrogenated, provided the
  topological information is provided.


Functions
---------

.. automodule:: mollib.hydrogens
    :members:
    :imported-members:
    :noindex:

Settings
--------

.. automodule:: mollib.hydrogens.settings
    :members:
    :noindex:


