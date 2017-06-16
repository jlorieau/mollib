Python Recipes
##############

Loading a Molecule
------------------

A molecule can be loaded directory from molecule objects. The following will
automatically download the file and cache it and read in the PDB for the
structure '2kxa'.

.. code:: python

    >>> from mollib import Molecule
    >>> mol = Molecule('2kxa')

This structure includes 10 models, and the first model is loaded by default.
Other models can be loaded.

.. code:: python

    >>> from mollib import Molecule
    >>> mol = Molecule('2kxa', model_id=3)

Loading Multiple Models
-----------------------

A list of models can be loaded using the :obj:`mollib.MoleculeReader` factory
directly.

.. code:: python

    >>> from mollib import MoleculeReader
    >>> mr = MoleculeReader()
    >>> all_models = mr.read('2kxa')
    >>> print(len(all_models))
    10

    >>> some_models = mr.read('2kxa', model_ids=[1, 5, 8])
    >>> for model in some_models:
    ...     print(model)
    Molecule (2kxa-1):    1 chains, 24 residues, 332 atoms.
    Molecule (2kxa-5):    1 chains, 24 residues, 332 atoms.
    Molecule (2kxa-8):    1 chains, 24 residues, 332 atoms.

Accessing Chains, Residues, Atoms and Coordinates
-------------------------------------------------

Molecules, chains, and residues are extented ``dict`` objects. Chains are
accessed by chain id (ex: ``'A'``), residues are accessed by residue number
(ex: ``3``), and atoms are accessed by atom name (ex: ``'CA'``).

.. code:: python

    >>> from mollib import Molecule
    >>> mol = Molecule('2kxa')  # loads the first model
    >>> chains = list(mol.chains)
    >>> print(chains)
    [A]

    >>> chainA = mol['A']
    >>> residues = list(chainA.residues)
    >>> print('{}, ...'.format(residues[:5]))
    [G1, L2, F3, G4, A5], ...

    >>> F3 = chainA[3]
    >>> atoms = list(F3.atoms)
    >>> print('{}, ...'.format(atoms[:5]))
    [A.F3.N, A.F3.CA, A.F3.C, A.F3.O, A.F3.CB], ...

    >>> CA = F3['CA']
    >>> print(round(CA.x, 1), round(CA.y, 1), round(CA.z, 1))
    (13.2, -2.7, 6.3)

    >>> mol['A'][3]['CA'] ==  CA
    True
