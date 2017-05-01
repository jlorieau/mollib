"""
The mollib is a simple Python toolset for manipulating molecular
structures in Python. A Molecule is a dict of Chains, which is a dict of
Residues, which is a dict of Atoms. The Molecule is constructed with helper
read and write functions in the Molecule class, and more sophisticated behavior
can be added to the Molecule, Chain, Residue and Atom classes by deriving them
and changing the chain_class, residue_class and atom_class of Molecule.

The :class:`Molecule` class accepts either a filename (in PDB format) or a
PDB id code, which it will download and cache. The new :obj:`molecule` object
is parsed in terms chains, residues and atoms.


>>> from mollib import Molecule
>>> mol=Molecule('2OED')
>>> print(mol)
Molecule (2OED):    1 chains, 56 residues, 862 atoms.


  .. note:: The :class:`Molecule` class will only read the first model for a
            structure that has many models.

Specific residues and their properties can be access. For example, the number
of atoms (:attr:`atom_size`) for residue #1 of chain 'A' is accessed using the
dict array notation.


>>> print("{} has {} atoms".format(mol['A'][1], mol['A'][1].atom_size))
M1 has 19 atoms


  .. note:: The number of atoms in a residue may be fewer than the actual number
            of atoms because some biomolecular structures do not include the
            positions of hydrogen atoms or atoms with large B-factors.


Likewise, the properties of specific atoms can be calculated. For example, the
mass and atomic position of :obj:`Atom` `CA` in residue 2:

>>> atom = mol['A'][2]['CA']
>>> print("Atom {} has a {:.2f} Da mass.".format(atom, atom.mass))
Atom A.Q2-CA has a 12.01 Da mass.
>>> print("{} is located at {}, {}, {} A".format(atom, *atom.pos))
A.Q2-CA is located at -3.086, -11.321, -1.361 A


Or properties of the entire molecule may be calculated.

>>> print("{} has a mass of {:.2f} Da".format(mol.name, mol.mass))
2OED has a mass of 6206.75 Da


The :class:`Molecule` class includes a number of methods to apply translations,
rotations and other operations on the molecule. The :attr:`center` method
centers the molecule.

>>> print("Center at {:.3f}, {:.3f}, {:.3f} A".format(*mol.center_of_mass))
Center at 0.133, -0.347, -0.002 A
>>> mol.center()
>>> center = map(abs, mol.center_of_mass)
>>> print("Center at {:.3f}, {:.3f}, {:.3f} A".format(*center))
Center at 0.000, 0.000, 0.000 A
"""

from . import settings
from .settings_manager import (register_settings, list_global_settings,
                               load_settings, import_settings)
from . import plugins
from .atom import Atom, sorted_atom_list
from .chain import Chain
from .geometry import *
from .molecule import Molecule
from .residue import Residue
from .utils import *

