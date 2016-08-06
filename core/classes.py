"""
MolLib for Python

MolLib is a simple Python toolset for manipulating molecular
structures in Python. A Molecule is a dict of Chains, which is a dict of
Residues, which is a dict of Atoms. The Molecule is constructed with helper
read and write functions in the Molecule class, and more sophisticated behavior
can be added to the Molecule, Chain, Residue and Atom classes by deriving them
and changing the chain_class, residue_class and atom_class of Molecule.

TODO: Add a Ensemble object that sits above Chain and below Molecule.
TODO: Replace dict inheritance with another object with dict-like accession.

>>> mol=Molecule('2OED')
>>> print(mol)
Molecule:    1 chains, 56 residues, 862 atoms.
>>> print("{} {:.0f}".format(mol['A'][1], mol['A'][1].atom_size))
M1 19
>>> print("{} {:.2f}".format(mol['A'][2]['CA'], mol['A'][2]['CA'].mass))
Q2-CA 12.01
>>> print("%.2f Da" % mol.mass)
6206.75 Da
>>> print("({:.3f}, {:.3f}, {:.3f})".format(*mol.center_of_mass))
(0.133, -0.347, -0.002)
>>> mol.rotate_zyz(0,90,0)
"""
# Author: Justin L Lorieau
# Copyright: 2016

import re
from itertools import chain as ichain
from math import cos, sin, sqrt, pi, atan2, acos
import numpy as np

from .utils import vector_length, calc_vector

# Imports for tests
import unittest
import doctest
from datetime import datetime

try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve

import os.path

# Utility Functions

# regex matching functions for string, floats and ints
re_str = re.compile(r'[a-zA-Z]')
re_float = re.compile(r'-?\d+\.\d*')
re_int = re.compile(r'-?\d+')


def convert(s):
    """Convert a string 's' into either an integer, float or string.
    Strips spaces.

    >>> value = convert('  -6.065')
    >>> print("{} {}".format(value, isinstance(value, float)))
    -6.065 True
    >>> value = convert('  3 ')
    >>> print("{} {}".format(value, isinstance(value, int)))
    3 True
    >>> value = convert(' 1232 ')
    >>> print("{} {}".format(value, isinstance(value, int)))
    1232 True
    >>> value = convert('HETATM  ')
    >>> print("{} {}".format(value, isinstance(value, str)))
    HETATM True
    >>> value = convert(' HZ3  ')
    >>> print("{} {}".format(value, isinstance(value, str)))
    HZ3 True
    """
    # If the string contains any letter, return as a string
    m = re_str.search(s)
    if m:
        return str(s).strip()

    # Try extracting float
    m = re_float.search(s)
    if m:
        return float(m.group())

    # Try extracting int
    m = re_int.search(s)
    if m:
        return int(m.group())

    # All else fails, try just returning the string
    return str(s).strip()


# MolLib Implementation

class PrimitiveMetaClass(type):
    """A Primitive Metaclass for properly assigning and collecting attributes,
    such as 'optional'"""
    def __new__(meta, classname, bases, classDict):

        # This code adds 'optional' tuples from parent classes to the
        # 'optional' attribute of this class.
        parent_optional = tuple(*[getattr(base, 'optional')
                                  for base in bases
                                  if hasattr(base, 'optional')])
        classDict['optional'] = classDict['optional'] + parent_optional \
            if 'optional' in classDict else parent_optional

        return type.__new__(meta, classname, bases, classDict)


class Primitive(object):
    "The base for all objects."

    __metaclass__ = PrimitiveMetaClass

    __slots__ = ()
    optional = ()

    def __init__(self, **kwargs):

        # Check that the required arguments have been specified
        req_kwargs = [kw for kw in self.__slots__ if kw not in self.optional]
        assert all(kw in kwargs for kw in req_kwargs), \
            "All of the following parameters are needed: {}"\
            .format(req_kwargs)

        # Assign the values
        [setattr(self, kw, value) for kw, value in kwargs.items()
         if kw in self.__slots__]

        super(Primitive, self).__init__()


class Atom(Primitive):
    """An atom in a residue.

    Parameters
    ----------
    number: int
        The atom number (order) in the molecule.
    name: str
        The atom's name. ex: 'HA'
    pos: :obj:`numpy.array`
        The atom's x/y/z position in the form of a numpy array.
    charge: float, optional
        The atom's charge
    element: str
        The atom's element. The element may be a different isotope. Suitable
        values are documented in atom_MW.
    residue: :obj:`Residue`, optional
        The Residue object to which this atom instance belongs to.
    chain: :obj:`Chain`, optional
        The Chain object to which this atom instance belongs to.
    molecule: :obj:`Molecule`, optional
        The Molecule object to which this atom instance belongs to.

    Attributes
    ----------
    atom_Mw: dict
        The molecular weight (value, in Da or g/mol) of each element type (key)
    options: tuple
        A list of fields that are optional
    """

    # These are the required field. 'pos' (position)is the coordinate position
    # of the atom, as a numpy array
    __slots__ = ('number', 'name', 'pos', 'charge', 'element',
                 'residue', 'chain', 'molecule')
    optional = ('charge', 'residue', 'chain', 'molecule')

    # Atom molecular weights. These must be labeled according the a str.title()
    # function. eg. ZN becomes Zn.
    atom_Mw = {'H': 1.01, 'C': 12.01, '13C': 13.00, 'N': 14.01, '15N': 15.00,
               'O': 16.00, 'Na': 22.99, 'Mg': 24.31, 'P': 30.97, 'S': 32.07,
               'Cl': 35.45, 'Zn': 65.38, }

    def __repr__(self):
        return u"{}-{}".format(self.residue, self.name) if self.residue else \
            u"{}".format(self.name)

    @property
    def mass(self):
        "The mass of this atom, depending on its element."
        return self.atom_Mw[self.element.title()]

    # Atom coordinate getters and setters
    @property
    def x(self):
        "The x-coordinate of this atom."
        return self.pos[0]

    @x.setter
    def x(self, value):
        self.pos[0] = value

    @property
    def y(self):
        "The y-coordinate of this atom."
        return self.pos[1]

    @y.setter
    def y(self, value):
        self.pos[1] = value

    @property
    def z(self):
        "The z-coordinate of this atom."
        return self.pos[2]

    @z.setter
    def z(self, value):
        self.pos[2] = value


class Residue(dict):
    """A residue in a chain.

    Parameters
    ----------
    name: str
        The residue's 3-letter name. ex: 'MET'
    letter: str
        The residue's 1-letter name. ex: 'M'
    number: int
        The residue's number in the sequence.

    Attributes
    ----------
    one_letter_codes: dict
        A dict for converting 3-letter residue names (key) to 1-letter residue
        names (value). The residue object will use 'X' if the 3-letter name
        is not found.
    last_residue: :obj:`Residue`
        The preceding residue object, in the sequence, from the molecule. This
        is a link in a doubly-linked list.
    next_residue: :obj:`Residue`
        The proceeding residue object, in the sequence, from the molecule. This
        is a link in a doubly-linked list.

        .. note:: The `last_residue` and `next_residue` are populated by the
                  Molecule object on creation--not the Residue object.
    """

    # These are linked-list pointers to the next and previous residues. These
    # attributes are populated during molecule creation
    last_residue = None
    next_residue = None

    one_letter_codes = {'ALA': 'A', 'GLY': 'G', 'SER': 'S', 'THR': 'T',
                        'MET': 'M', 'CYS': 'C', 'ILE': 'I', 'LEU': 'L',
                        'VAL': 'V', 'PHE': 'F', 'TYR': 'Y', 'TRP': 'W',
                        'ASN': 'N', 'GLN': 'Q', 'ASP': 'D', 'GLU': 'E',
                        'HIS': 'H', 'PRO': 'P', 'ARG': 'R', 'LYS': 'K'}

    def __init__(self, name, number, *args, **kwargs):
        name = str(name).upper()

        self.name = name  # full name, MET
        self.letter = self.one_letter_codes.get(self.name, 'X')  # letter code
        self.number = number
        super(Residue, self).__init__(*args, **kwargs)

    def __repr__(self):
        return u"{}{}".format(self.letter, self.number)

    @property
    def atoms(self):
        "An iterator over all atoms in this residue, sorted by atom number"
        return (atom for atom in sorted(self.values(), key=lambda a: a.number))

    @property
    def atom_size(self):
        "The number of atoms in this residue."
        return len(list(self.atoms))

    @property
    def mass(self):
        """The mass of all atoms in this residue.

        >>> mol = Molecule('2KXA')
        >>> print("{:.2f}".format(mol['A'][3].mass)) # Phe-3 mass
        147.19
        """
        return sum(a.mass for a in self.atoms)

    @property
    def ramachandran_angles(self):
        """The backbone Ramachandran angles for this residue.

        Returns
        -------
        angles: list
            A pair of angles (floats) representing the phi and psi angles in
            degrees.

        >>> mol = Molecule('2KXA')
        >>> print("{:.1f} {:.1f}".format(*mol['A'][3].ramachandran_angles))
        -65.2 -48.3
        >>> print("{} {:.1f}".format(*mol['A'][1].ramachandran_angles))
        None -179.9
        """
        # Retrieve the appropriate atoms to calculate the ramachandran angles
        c_prev = (self.last_residue['C']
                  if self.last_residue is not None else None)
        n = self['N']
        ca = self['CA']
        c = self['C']
        n_next = (self.next_residue['N']
                  if self.next_residue is not None else None)

        angles = [None, None]  # Initial phi/psi angles are set to None
        # Iterate over phi (angles[0]) and phi (angles[1]) to calculate the
        # dihedrals
        for item, a, b, c, d in [(0, c_prev, n, ca, c),   # phi
                                 (1, n, ca, c, n_next)]:  # psi
            # The atom has to exist to calculate a dihedral
            if a is None or b is None or c is None or d is None:
                continue

            # Calculate the normalized vectors.
            ab = calc_vector(a, b, normalize=True)
            bc = calc_vector(b, c, normalize=True)
            cd = calc_vector(c, d, normalize=True)

            # The dihedral is the angle between the plans a-b-c and b-c-d
            # The angle between these plans can be calculated from their
            # normals (cross products)
            n1 = np.cross(ab, bc)
            n2 = np.cross(bc, cd)

            # The angle between n1 and n2 can be calculated with acos. However
            # the followin atan2 relationship returns a number between 0 and
            # 2pi
            m1 = np.cross(n1, bc)
            x = np.dot(n1, n2)
            y = np.dot(m1, n2)
            angle = atan2(y, x) * 180. / np.pi

            angles[item] = angle

        return angles


class Chain(dict):
    """A chain in a molecule.

    Parameters
    ----------
    id: str
        The chain's id. ex: 'A'
    """

    def __init__(self, id, *args, **kwargs):
        self.id = id
        super(Chain, self).__init__(*args, **kwargs)

    def __repr__(self):
        return u"{}".format(self.id)

    @property
    def residues(self):
        """An iterator over all residues in this chain,
        sorted by residue number."""
        return (residue for residue in
                sorted(self.values(), key=lambda a: a.number))

    @property
    def residue_size(self):
        "The number of residues in this chain."
        return len(list(self.residues))

    @property
    def atoms(self):
        """An iterator over all atoms in this chain,
        sorted by atom number."""
        return (a for a in sorted(ichain(*[r.values() for r in self.values()]),
                                  key=lambda a: a.number))

    @property
    def atom_size(self):
        "The number of atoms in this chain."
        return len(list(self.atoms))

    @property
    def mass(self):
        """The mass of all the atoms in this chain.

        >>> mol = Molecule('3C9J')
        >>> print("{:.2f}".format(mol['A'].mass)) # Chain A mass
        2507.61
        >>> print("{:.2f}".format(mol.mass)) # Total molecule mass
        10164.55
        """
        return sum(a.mass for a in self.atoms)


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

        .. note:: A molecule object only reads the first model of a PDB file with
                  multiple models
    """

    # TODO: add translate method
    # TODO: Add atom notes
    # TODO: Add atom isotopes

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

        >>> mol=Molecule('2KXA')
        >>> print([r.number for r in mol.residues])
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, \
20, 21, 22, 23, 24]
        """
        return (r for r in sorted(ichain(*[r.values() for r in self.values()]),
                                  key=lambda r: r.number))

    @property
    def residue_size(self):
        "The number of residues in this molecule."
        return len(list(self.residues))

    @property
    def atoms(self):
        "An iterator over all atoms in this molecule, sorted by atom number"
        return (a for a in sorted(ichain(*[r.values() for r in
                                           ichain(*[c.values() for c in
                                                    self.values()])]),
                                  key=lambda a: a.number))

    @property
    def atom_size(self):
        "The number of atoms in this molecule."
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
        "Rotates a molecule by the Euler z-y-z angles in degrees."

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

        atom = self.atom_class(**kwargs)
        atom.residue = residue
        atom.chain = residue.chain
        atom.molecule = self
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
        "Create a doubly linked list of all residues."
        # Create the residue linked lists
        last_residue = None
        next_residue = None
        for residue in self.residues:
            # Special first residue
            if last_residue is None:
                last_residue = residue
                continue
            last_residue.next_residue = residue
            residue.last_residue = last_residue
            last_residue = residue

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
        "Reads in data from a PDB file."

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
        url = 'http://ftp.rcsb.org/download/{}.pdb'.format(pdb_code)
        path = os.path.join('/tmp', pdb_code) + '.pdb'

        if not os.path.isfile(path) or not load_cached:
            urlretrieve(url, path)
        self.read_pdb(path)

    def read_stream(self, stream):
        "Reads in data from a stream."

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
                match = pdb_line.match(line)
                if match:
                    yield match

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
            id = groupdict['chain']

            # If this is a HETATM, then append a '*' to the chain name so that
            # it doesn't overwrite protein chains.
            if groupdict['type'] == 'HETATM':
                id += '*'

            if id not in self:
                chain = self.chain_class(id=id)
                chain.molecule = self
                self[id] = chain

            chain = self[id]

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
            atom_dict = {k: v for k, v in groupdict.items()
                         if k in Atom.__slots__}
            atom = self.atom_class(**atom_dict)
            atom.residue = residue
            atom.chain = chain
            atom.molecule = self
            residue[name] = atom

        self.link_residues()

    # Molecular Properties

    @property
    def mass(self):
        """ The mass of the molecule.

        Examples
        --------
        >>> mol = Molecule('2KXA')
        >>> print("{:.2f}".format(mol.mass))
        2445.07
        """
        return sum(a.mass for a in self.atoms)

    @property
    def center_of_mass(self):
        """The center-of-mass position of the molecule.

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


# TESTS #

class TestMolLib(unittest.TestCase):

    performance_tests = False

    def test_large_molecule(self):
        "Tests the parsing and performance of a very large protein complex."
        import string

        if self.performance_tests:
            # 5.5s - default atom generator
            # 6.0s - swithching to atom generator 2
            # 3.6s - after switching convert function to regexes. (2.75MB/s)
            id = '3H0G'  # RNA Polymerase II from Schizosaccharomyces pombe
            start = datetime.now()
            Molecule(id)
            stop = datetime.now()
            time = (stop - start).total_seconds()
            print("Loaded {id} in {time} seconds".format(id=id, time=time))

        mol = Molecule('3H0G')

        # Test that all of the chains were read in correctly
        chains = list(string.ascii_uppercase)[:24]  # chains A-X
        chains += ['A*', 'B*', 'C*', 'I*', 'J*', 'L*', 'M*', 'N*', 'O*', 'U*',
                   'V*', 'X*']
        self.assertEqual([c.id for c in mol.chains], sorted(chains))
        self.assertEqual(mol.chain_size, len(chains))

        # Test the molecular mass of each chain
        for chain in mol.chains:
            self.assertGreater(chain.mass, 0.)

        self.assertAlmostEqual(mol.mass, 833388.28, 2)

    def test_multiple_models(self):
        """Tests reading PDB files with multiple models. Only the first model
        should be read in."""
        mol = Molecule('2KXA')  # 20 models

        # These are the coordinates for this atom of the first model
        self.assertEqual(mol['A'][3]['N'].pos[0], 13.766)
        self.assertEqual(mol['A'][3]['N'].pos[1], -3.965)
        self.assertEqual(mol['A'][3]['N'].pos[2], 5.893)

    def test_residue_ordering(self):
        """Tests the linked lists of residues."""
        mol = Molecule('2KXA')

        last_residues = [r.last_residue.number
                         if r.last_residue is not None else None
                         for r in mol.residues]
        self.assertEqual(last_residues,
                         [None, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                          14, 15, 16, 17, 18, 19, 20, 21, 22, 23])
