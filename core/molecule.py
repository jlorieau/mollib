"""
MolLib Overview
===============
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
Molecule:    1 chains, 56 residues, 862 atoms.


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
Atom Q2-CA has a 12.01 Da mass.
>>> print("{} is located at {}, {}, {} A".format(atom, *atom.pos))
Q2-CA is located at -3.086, -11.321, -1.361 A


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

# Author: Justin L Lorieau
# Copyright: 2016
# TODO: Add a Ensemble object that sits above Molecule.
# TODO: Add atom annotations, like hbond_donor and hydridization (sp2, sp3)
#       hbond_donor = True, hbond_acceptor = True, hybridization = 'sp2'

import re
import weakref
from itertools import chain as ichain
from math import cos, sin, pi

import numpy as np

from .utils import convert
from .atom import Atom
from .residue import Residue
from .chain import Chain
from . import settings

try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve

import os.path


class Molecule(dict):
    """A class for molecular structures.

    Parameters
    ----------
    name: str
        The name of this molecule. This is usually the PDB code of the molecule.

    Attributes
    ----------
    chain_class: :class:`Chain`
        The `Chain` class to use for generating chain objects for this
        molecule.
    residue_class: :class:`Residue`
        The `Residue` class to use for generating residue objects for this
        molecule.
    atom_class: :class:`Atom`
        The `Atom` class to use for generating atom objects for this molecule.

        .. note:: A molecule object only reads the first model of a PDB file
                  with multiple models
    """

    # TODO: add translate method
    # TODO: Add atom notes
    # TODO: Add atom isotopes
    # TODO: Add function to read pdb header

    # The following class-level attributes are used to customize the base or
    # derived Chain, Residue and Atom classes used in the molecule
    chain_class = Chain
    residue_class = Residue
    atom_class = Atom

    def __init__(self, identifier, *args, **kwargs):
        """Molecule constructor that accepts an identifier.

        Parameters
        ----------
        identifier: str
            An molecule identifier that is either a filename (PDB format),
            or PDB code (ex: '2KXA').
        """
        self.name = identifier
        self.read_identifier(identifier)
        super(Molecule, self).__init__(*args, **kwargs)

    def __repr__(self):
        return (u"Molecule:"
                "    {} chains, {} residues, {} atoms."
                .format(self.chain_size, self.residue_size, self.atom_size))

    # Basic Accessors and Mutators

    @property
    def chain_size(self):
        """The number of chains in this molecule.

        >>> mol=Molecule('1HTM') # 6 subunits, 4 types of HOH
        >>> print(mol.chain_size)
        10
        """
        return len(self)

    @property
    def chains(self):
        """An iterator over all chains in this molecule, sorted by chain id.

        >>> mol=Molecule('3C9J')
        >>> print([c.id for c in mol.chains])
        ['A', 'B', 'B*', 'C', 'D']
        """
        return (c for k, c in sorted(self.items()))

    @property
    def residues(self):
        """An iterator over all residues in this molecule, sorted by residue
        number.

        >>> mol=Molecule('2N65')
        >>> l = [r.number for r in mol.residues]
        >>> print(l[16:])
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        >>> print(l[:16])
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        """
        return (r for r in sorted(ichain(*[r.values() for r in self.values()]),
                                  key=lambda r: (r.chain.id, r.number)))

    @property
    def residue_size(self):
        """The number of residues in this molecule."""
        return len(list(self.residues))

    @property
    def atoms(self):
        """An iterator over all atoms in this molecule, sorted by atom
        number."""
        return (a for a in sorted(ichain(*[r.values() for r in
                                           ichain(*[c.values() for c in
                                                    self.values()])]),
                                  key=lambda a: a.number))

    @property
    def atom_size(self):
        """The number of atoms in this molecule."""
        return len(list(self.atoms))

    # Mutator Functions

    def center(self):
        """Centers a molecule about its center_of_mass.

        >>> mol = Molecule('2KXA')
        >>> print("{:.3f} {:.3f} {:.3f}".format(*mol.center_of_mass))
        16.938 -0.058 0.125
        >>> mol.center()
        >>> # Center molecule. Map to absolute value to avoid values of -0.000
        >>> print("{:.3f} {:.3f} {:.3f}".format(*map(abs, mol.center_of_mass)))
        0.000 0.000 0.000
        """
        com = self.center_of_mass
        for atom in self.atoms:
            atom.pos[0] -= com[0]
            atom.pos[1] -= com[1]
            atom.pos[2] -= com[2]

    def rotate_zyz(self, alpha, beta, gamma):
        """Rotates a molecule by the Euler z-y-z angles in degrees."""

        # Trig.
        sin_a = sin(alpha * pi / 180.)
        sin_b = sin(beta * pi / 180.)
        sin_g = sin(gamma * pi / 180.)

        cos_a = cos(alpha * pi / 180.)
        cos_b = cos(beta * pi / 180.)
        cos_g = cos(gamma * pi / 180.)

        m = np.matrix([[-sin_a * sin_g + cos_a * cos_b * cos_g,
                        cos_a * sin_g + sin_a * cos_b * cos_g,
                        -sin_b * cos_g],
                       [-sin_a * cos_g - cos_a * cos_b * sin_g,
                        cos_a * cos_g - sin_a * cos_b * sin_g,
                        sin_b * sin_g],
                       [cos_a * sin_b,
                        sin_a * sin_b,
                        cos_b]])
        for atom in self.atoms:
            v = np.matrix([atom.pos[0], atom.pos[1], atom.pos[2]]).T
            v_new = np.dot(m, v)
            v_new = v_new.tolist()
            atom.pos = np.array((v_new[0][0], v_new[1][0], v_new[2][0]))

        return None

    def add_atom(self, name, pos, charge, element, residue, **kwargs):
        """Adds an atom to the molecule.

        Parameters
        ----------
        name: str
            The name of the atom. ex: 'HA'
        pos: :obj:`numpy.array`
            An array containing the x, y and z coordinates of the atom.
        charge: float
            The charge of the atom.
        element: str
            The type of element for the atom. ex: 'H'
        residue: :obj:`Residue`
            The residue object to add the item to. This molecule should already
            contain the residue object.
        kwargs: dict, optional
            Additional parameters to pass to the Atom constructor.

        Examples
        --------
        >>> mol = Molecule('2KXA')
        >>> print ('{:.2f}, {:.2f}'.format(mol.mass, mol['A'][3].mass))
        2445.07, 147.19
        >>> mol.add_atom('C3', (0.0, 0.0, 0.0), 0., 'C', mol['A'][3])
        >>> print ('{:.2f}, {:.2f}'.format(mol.mass, mol['A'][3].mass))
        2457.08, 159.20
        """
        kwargs['number'] = -1
        kwargs['name'] = name
        kwargs['pos'] = np.array(pos)
        kwargs['charge'] = charge
        kwargs['element'] = element
        kwargs['residue'] = residue
        kwargs['chain'] = residue.chain
        kwargs['molecule'] = self

        atom = self.atom_class(**kwargs)
        residue[name] = atom

    def del_atom(self, atom):
        """Deletes the specified :obj:`atom` from the molecule.

        Examples
        --------
        >>> mol = Molecule('2KXA')
        >>> print ('{:.2f}, {:.2f}'.format(mol.mass, mol['A'][3].mass))
        2445.07, 147.19
        >>> mol.del_atom(mol['A'][3]['N'])
        >>> print ('{:.2f}, {:.2f}'.format(mol.mass, mol['A'][3].mass))
        2431.06, 133.18


        .. note:: Since atom objects are only referenced in residues, this
                  should delete the only atom reference in the molecule.
                  However, if the atom was referenced elsewhere, it won't get
                  deleted until that reference is deleted. An implementation
                  that returns weakref.proxy links from the molecule may be
                  preferable.
        """
        del atom.residue[atom.name]

    def strip_atoms(self, element):
        """Deletes all atoms with the given element (name, i.e. 'H')

        Examples
        --------
        >>> mol = Molecule('2KXA')
        >>> print('{:.2f}'.format(mol.mass))
        2445.07
        >>> mol.strip_atoms('H')
        >>> print('{:.2f}'.format(mol.mass))
        2285.49
        """
        for atom in self.atoms:
            if atom.element == element:
                del atom.residue[atom.name]

    def substitute_element(self, element_from, element_to):
        """Substitutes all atoms of one element to another.

        Examples
        --------
        >>> mol = Molecule('2KXA')
        >>> print('{:.2f}'.format(mol.mass))
        2445.07
        >>> mol.substitute_element('C', '13C')
        >>> print('{:.2f}'.format(mol.mass))
        2559.91
        """
        for atom in self.atoms:
            if atom.element == element_from:
                atom.element = element_to

    def link_residues(self):
        """Create a doubly linked list of all residues and annotate the
        first and last residues."""
        # Create the residue linked lists
        # Note: the prev_residue, next_residue, first and last attributes
        # of the residue object are already set by the class
        prev_residue = None
        residue = None
        for residue in self.residues:
            # Treat the first residue in the chain as special
            if prev_residue is None:
                residue.first = True

                # Set the double-linked list
                residue.prev_residue = None

            # Otherwise this is not the first residue
            else:
                # Treat the first residue in the *next* chain as special
                if prev_residue.chain.id != residue.chain.id:
                    residue.first = True
                    prev_residue.last = True

                    # Set the double-linked list
                    prev_residue.next_residue = None
                    residue.prev_residue = None
                else:
                    # Set the double-linked list
                    prev_residue.next_residue = weakref.proxy(prev_residue)
                    residue.prev_residue = weakref.proxy(residue)

            # Prepare for the next iteration
            prev_residue = residue

        # Treat the very last residue as the last residue
        if residue is not None:
            residue.last = True
            residue.next_residue = None

    @property
    def pH(self):
        return getattr(self, '_pH', settings.default_pH)

    @pH.setter
    def pH(self, value):
        self._pH = value

    # Read and Write Methods

    def write_pdb(self, filename):
        """Write data to a PDB file.

        Parameters
        ----------
        filename : str
            The path and filename to write the file to.
        """

        with open(filename, 'w') as f:

            # Populate 'ATOM' lines

            atom_line = ("{line_type:6}"
                         "{atom_num:>5} "
                         "{atom_name:<4}"
                         "{alt_loc:1}"
                         "{res_name:3} "
                         "{chain:1}"
                         "{res_number:4}"
                         "{icode:1}   "
                         "{x:8.3f}"
                         "{y:8.3f}"
                         "{z:8.3f}"
                         "{occupancy:6.2f}"
                         "{B_factor:6.2f}"
                         "          "
                         "{element:>2}"
                         "{charge:2}"
                         "\n")

            for count, atom in enumerate(self.atoms, 1):
                atom_parms = {'line_type': 'ATOM',
                              'atom_num': count,
                              'atom_name': atom.name,
                              'alt_loc': '',
                              'res_name': atom.residue.name,
                              'chain': atom.chain,
                              'res_number': atom.residue.number,
                              'icode': '',
                              'x': atom.pos[0],
                              'y': atom.pos[1],
                              'z': atom.pos[2],
                              'occupancy': 1,
                              'B_factor': 0,
                              'element': atom.element,
                              'charge': atom.charge}
                f.write(atom_line.format(**atom_parms))

    def read_identifier(self, identifier):
        """Reads in structure based on an identifier

        Parameters
        ----------
        identifier : str
            The `identifier` is either a filename, path, or 4-alphanumeric
            PDB code.
        """
        # Check to see if identifier is a filename or path
        if os.path.isfile(identifier):
            self.read_pdb(identifier)
        else:
            self.fetch_pdb(identifier)

    def read_pdb(self, filename):
        """Reads in data from a PDB file."""

        with open(filename) as f:
            self.read_stream(f)

    def fetch_pdb(self, pdb_code, load_cached=True):
        """Download/fetch a PDB file online.

        Parameters
        ----------
        pdb_code : str
            The 4 alphanumeric character PDB code to load.
        load_cached : bool, optional
            If a cached version is available, use that instead of downloading
            the file.

        Examples
        --------
        >>> mol = Molecule('2KXA')
        >>> mol.fetch_pdb(pdb_code='2KXA', load_cached=False) # force download
        >>> print(mol)
        Molecule:    1 chains, 24 residues, 332 atoms.
        """
        # TODO: fetch gzipped files instead of raw PDB files.
        url = 'http://ftp.rcsb.org/download/{}.pdb'.format(pdb_code)
        path = os.path.join('/tmp', pdb_code) + '.pdb'

        if not os.path.isfile(path) or not load_cached:
            urlretrieve(url, path)
        self.read_pdb(path)

    def read_stream(self, stream):
        """Reads in data from a stream.
        """
        # TODO: This may be much faster to implement in Cython with fixed str.
        self.clear()

        pdb_line = re.compile((r"(?P<type>ATOM  |HETATM)"
                               "(?P<number>[\s\d]{5}) "
                               "(?P<name>[\s\w]{4})"
                               "(?P<alt_loc>[\w\s])"
                               "(?P<residue_name>[\w\s]{3}) "
                               "(?P<chain>[\s\w]{1})"
                               "(?P<residue_number>[\s\w]{4})"
                               "(?P<icode>[\w\s])   "
                               "(?P<x>[\d\s\.\-]{8})"
                               "(?P<y>[\d\s\.\-]{8})"
                               "(?P<z>[\d\s\.\-]{8})"
                               "(?P<occupancy>[\d\s\.\-]{6})"
                               "(?P<B_factor>[\d\s\.\-]{6})          "
                               "(?P<element>[\s\w]{2})"
                               "(?P<charge>[\d\s\.\-]{2})?"))

        # Find the ATOM/HETATM lines and pull out the necessary data

        # Generator implementation 1: This generator function reads only the
        # the first model. It's a little slower than implementation 2 (about
        # 10%). It takes 6.0s to read 3H0G.
        def generator():
            for line in stream.readlines():
                # Read only the first model. Skip the rest
                if line[0:6] == 'ENDMDL':
                    raise StopIteration
                m = pdb_line.match(line)
                if m:
                    yield m

        atom_generator = generator()

        # Generator implementation 2: This generator is a little faster than
        # generator 1, but it reads all of the models, saving only the last
        # one. It takes 5.5 seconds on 5H0G
        # atom_generator = filter(None, map(pdb_line.match,
        #                                  stream.readlines()))

        # Retrieve a set from the match objects
        for match in atom_generator:
            groupdict = {field_name: convert(field_value)
                         for field_name, field_value
                         in match.groupdict().items()}

            # create Chain, if it doesn't already exist
            identifier = groupdict['chain']

            # If this is a HETATM, then append a '*' to the chain name so that
            # it doesn't overwrite protein chains.
            if groupdict['type'] == 'HETATM':
                identifier += '*'

            # Create a new chain, if it doesn't already exist
            if identifier not in self:
                chain = self.chain_class(identifier=identifier)
                chain.molecule = self
                self[identifier] = chain

            chain = self[identifier]

            # create Residue, if it doesn't already exist
            number, name = (groupdict[i] for i in ('residue_number',
                                                   'residue_name'))
            if number not in chain:
                try:
                    residue = self.residue_class(number=number, name=name)

                    residue.chain = chain
                    residue.molecule = self
                    chain[number] = residue
                except KeyError:
                    continue
            residue = chain[number]

            # create the Atom. The following code overwrites atoms duplicate
            # in atom name
            name = groupdict['name']

            # Reformat the x/y/z coordinates to a numpy array
            groupdict['pos'] = np.array((groupdict.pop('x'),
                                         groupdict.pop('y'),
                                         groupdict.pop('z')))

            # Populate the new atom parameters and create the new Atom object.
            atom_dict = {k: v for k, v in groupdict.items()
                         if k in Atom.__slots__}
            atom_dict['residue'] = residue
            atom_dict['chain'] = chain
            atom_dict['molecule'] = self

            atom = self.atom_class(**atom_dict)
            residue[name] = atom

        self.link_residues()

    # Molecular Properties

    @property
    def mass(self):
        """Mass of the molecule.

        Examples
        --------
        >>> mol = Molecule('2KXA')
        >>> print("{:.2f}".format(mol.mass))
        2445.07
        """
        return sum(a.mass for a in self.atoms)

    @property
    def center_of_mass(self):
        """Center-of-mass position of the molecule.

        Returns
        -------
        pos : :obj:`numpy.array`
            The x, y and z coordinates.

        Examples
        --------
        >>> mol = Molecule('2KXA')
        >>> print("{:.3f} {:.3f} {:.3f}".format(*mol.center_of_mass))
        16.938 -0.058 0.125
        """

        pos = np.array((0., 0., 0.))
        m_total = 0.
        for atom in self.atoms:
            mass = atom.mass
            m_total += mass
            pos += atom.pos*mass
        return pos / m_total
