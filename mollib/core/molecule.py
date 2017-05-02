"""
The molecule object.
"""

# Author: Justin L Lorieau
# Copyright: 2016
# TODO: Separate pdb read and writer to separate file
# TODO: Make PDB read parallel
# TODO: Improve PDB writer
# TODO: Add a Ensemble object that sits above Molecule.
# TODO: Add atom annotations, like hbond_donor and hydridization (sp2, sp3)
#       hbond_donor = True, hbond_acceptor = True, hybridization = 'sp2'

import re
import weakref
import logging
import gzip
from itertools import chain as ichain
from itertools import count
from collections import OrderedDict
from math import cos, sin, pi
import os.path

import numpy as np

from .atom import Atom
from .residue import Residue
from .chain import Chain
from .topology import topology
from .utils import grouper
from mollib.utils.net import get_or_fetch
from . import settings


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
    connections: list
        A list of lists of atom numbers for atoms that are connected to each
        other. Note that each list item respect the field order and position
        of the official PDB format.

    _parameters: dict
        Stores the molecule's parameters
    cache: dict
        Stored cached calculations for the molecule. These are cleared by
        the clear_cache method.

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


    def __new__(cls, *args, **kwargs):
        "Keep track of class instances"
        instance = dict.__new__(cls, *args, **kwargs)
        if "_instances" not in cls.__dict__:
            cls._instances = []
        ref = weakref.ref(instance)
        if ref not in cls._instances and ref is not None:
            cls._instances.append(ref)
        return instance

    def __init__(self, identifier, *args, **kwargs):
        """Molecule constructor that accepts an identifier.

        Parameters
        ----------
        identifier: str
            An molecule identifier that is either a filename (PDB format),
            or PDB code (ex: '2KXA').
        """
        # The following strips path and extensition information from the
        # identifier to make an easily readable name.
        name = os.path.split(identifier)[-1]
        name = os.path.splitext(name)[0]

        self.name = name
        self.identifier = identifier
        self.connections = []
        self._parameters = {}
        self.cache = {}

        # Read in the data
        self.read_identifier(identifier)
        super(Molecule, self).__init__(*args, **kwargs)

    def __repr__(self):
        return ("Molecule ({}):".format(self.name) +
                "    {} chains, {} residues, {} atoms."
                .format(self.chain_size, self.residue_size, self.atom_size))

    # Class properties

    @classmethod
    def get_weakref(cls, name, instance_number=0):
        """Return a weakreference to an *existing* :obj:`molecule` instance
        with the given name.

        Parameters
        ----------
        name: str
            The name of the molecule. ex: '2KXA'

        Returns
        -------
        weakref.ref
            A weakreference to the molecule instance.
        None
            If the molecule was not found.
        """
        if (hasattr(cls, '_instances')):
            refs = [i for i in cls._instances
                    if i() is not None and i().name == name]
            if len(refs) > 0:
                return refs[0]
        return None

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

        Examples
        --------
        >>> mol=Molecule('2N65')
        >>> l = [r.number for r in mol.residues]
        >>> print(l[16:])
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        """
        return (r for r in sorted(ichain(*[r.values() for r in self.values()]),
                                  key=lambda r: (r.chain.id, r.number)))

    @property
    def residues_reversed(self):
        """An iterator over all residues in this molecule, sorted in *reversed*
        residue number.

        Examples
        --------
        >>> mol=Molecule('2N65')
        >>> l = [r.number for r in mol.residues_reversed]
        >>> print(l[16:])
        [16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        """
        return (r for r in sorted(ichain(*[r.values() for r in self.values()]),
                                  key=lambda r: (r.chain.id, r.number),
                                  reverse=True))

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

    _selector = re.compile(r'\s*'
                           '((?P<chains>[A-Z]+\*?:?[A-Z]*\*?)\.)?'
                           '(?P<res_nums>\d+:?\d*)'
                           '\.?'
                           '(?P<atom_name>[A-Z]+\d*)\s*')

    def get_atoms(self, *selectors):
        """Return atoms matches the given locator names.

        Parameters
        ----------
        selectors:
            A series of selector strings that match self._selector. The atom
            locators use the following conventions:

            1. (residue number)-(atom name). ex: 31-CA
            2. (chain id)-(residue number)-(atom name). ex: A.31-CA

            .. note::  Residue number of chain id ranges are allowed using a
                       ':'. ex: 18:24-CA or A:D.13-CB.

            If no chain id is specified, the chain 'A' is used.

        Returns
        -------
        list
            A list of :obj:`mollib.Atom` objects.

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2MUV')
        >>> mol.get_atoms('A.22.CA', '26.CB')
        [A.S22.CA, A.L26.CB]
        >>> mol.get_atoms('A.335-CA')  # doesn't exist
        []
        >>> mol.get_atoms('A.22:26.CA')
        [A.S22.CA, A.S23.CA, A.D24.CA, A.P25.CA, A.L26.CA]
        >>> mol.get_atoms('A:D.22.CA')
        [A.S22.CA, B.S22.CA, C.S22.CA, D.S22.CA]
        """
        # TODO: Add first and last residue selectors. (f:l)
        atoms = []

        # Find the atoms given by the names
        for selector in selectors:
            # Find the relevant atoms using the locator
            match = self._selector.match(selector)

            if match is None:
                msg = ("The selector '{}' could not be interpreted.")
                logging.error(msg.format(selector))
                continue

            # Parse the selector
            seldict = match.groupdict()

            # Convert A:D into ('A', 'B', 'C', 'D')
            chain_ids = seldict['chains']
            chain_ids = 'A' if chain_ids is None else chain_ids
            chain_ids = chain_ids.split(':')
            if len(chain_ids) == 2:
                id1, id2 = chain_ids
                chain_ids = [chr(i) for i in range(ord(id1), ord(id2) + 1)]

            # Convert 23:26 into (23, 24, 25, 26) and convert into integers
            residue_numbers = seldict['res_nums']
            residue_numbers = residue_numbers.split(':')
            if len(residue_numbers) == 2:
                num1, num2 = int(residue_numbers[0]), int(residue_numbers[1])
                residue_numbers = [int(i) for i in range(num1, num2 + 1)]
            else:
                residue_numbers = [int(r) for r in residue_numbers]

            atom_name = seldict['atom_name']

            # Get the relevant chains, residues and (hopefully) atoms
            for chain_id in chain_ids:
                chain = self.get(chain_id, None)
                if chain is None:
                    continue

                for residue_number in residue_numbers:
                    residue = chain.get(residue_number, None)
                    if residue is None:
                        continue

                    atom = residue.get(atom_name, None)

                    if atom is None:
                        continue

                    atoms.append(atom)

        if not atoms:
            msg = ("The selector '{}' could not find any atoms.")
            logging.error(msg.format(selector))
        return atoms

    # Getting and setting parameters
    def get_parameter(self, category, name):
        """Get the parameter for the given property category and name.

        Parameters
        ----------
        category: str
            The category of the parameter. ex: 'hydrogens'
        name: str
            The name of the parameter. ex: 'A.Q61.CA'

        Returns
        -------
        value
            The property value for the name parameter in the given category.
        None
            If no value was found.

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2KXA')
        >>> mol.set_parameter('Biophysical Characteristics', 'Rg', 10.5)
        >>> mol.get_parameter('Biophysical Characteristics', 'Rg')
        10.5
        """
        category_dict = self._parameters.setdefault(category, {})
        return category_dict.get(name, None)

    def set_parameter(self, category, name, value):
        """Set the parameter for the given property category and name.

        Parameters
        ----------
        category: str
            The category of the parameter.
        name: str
            The name of the parameter. ex: 'A.Q61.CA'
        value:
            The value of the parameter


        .. note:: Common categories include

                  - 'Structural Features': These are a listing of parameters
                    that highlight structural features within the molecule.
                    They're expected to be invariant to whole-body translations
                    and rotations, but not changes in internal structure.
        """
        category_dict = self._parameters.setdefault(category, {})
        category_dict[name] = value

    # Mutator Functions

    def center(self):
        """Centers a molecule about its center_of_mass.


        .. note:: This function invalidates cache objects with the attribute
                  or key 'preserve_cache_wb_translation'

        Examples
        --------
        >>> mol = Molecule('2KXA')
        >>> print("{:.3f} {:.3f} {:.3f}".format(*mol.center_of_mass))
        16.938 -0.058 0.125
        >>> mol.center()
        >>> # Center molecule. Map to absolute value to avoid values of -0.000
        >>> print("{:.3f} {:.3f} {:.3f}".format(*map(abs, mol.center_of_mass)))
        0.000 0.000 0.000
        """
        # Invalidate the caches
        self.clear_cache(scope='wb_translation')

        com = self.center_of_mass
        for atom in self.atoms:
            atom.pos[0] -= com[0]
            atom.pos[1] -= com[1]
            atom.pos[2] -= com[2]

    def rotate_zyz(self, alpha, beta, gamma):
        """Rotates a molecule by the Euler z-y-z angles in degrees.

        .. note:: This function invalidates the caches without
                  'preserve_cache_wb_rotation'
        """
        # TODO: Use the rotation function in the utils
        # Invalidate the caches
        self.clear_cache(scope='wb_rotation')

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

    def add_atom(self, name, pos, element, residue, bonded_atoms = None,
                 clear_cache=True, **kwargs):
        """Adds an atom to the molecule.

        Parameters
        ----------
        name: str
            The name of the atom. ex: 'HA'
        pos: :obj:`numpy.array`
            An array containing the x, y and z coordinates of the atom.
        element: str
            The type of element for the atom. ex: 'H'
        residue: :obj:`Residue`
            The residue object to add the item to. This molecule should already
            contain the residue object.
        bonded_atoms: list
            If specified, these atoms will be added to this atom's topology
            (and vice-versa)
        clear_cache: bool, optional
            If True, this function will clear the cache for objects. See note
            below
        kwargs: dict, optional
            Additional parameters to pass to the Atom constructor.


        .. note:: This function invalidates cache objects without the attribute
                  or method 'preserve_cache_add_atoms'

        Examples
        --------
        >>> mol = Molecule('2KXA')
        >>> print ('{:.2f}, {:.2f}'.format(mol.mass, mol['A'][3].mass))
        2445.07, 147.19
        >>> mol.add_atom('C3', (0.0, 0.0, 0.0), 'C', mol['A'][3])
        >>> print ('{:.2f}, {:.2f}'.format(mol.mass, mol['A'][3].mass))
        2457.08, 159.20
        """
        # Invalidate the wb_rotation cache
        if clear_cache:
            self.clear_cache(scope='add_atoms')

        # TODO: add test.
        if 'number' not in kwargs:
            kwargs['number'] = -1
        # kwargs['number'] = -1
        kwargs['name'] = name
        kwargs['pos'] = np.array(pos)
        kwargs['element'] = element
        kwargs['residue'] = residue
        kwargs['chain'] = residue.chain
        kwargs['molecule'] = self

        atom = self.atom_class(**kwargs)
        residue[name] = atom

        # Add this atom to the topologies
        if hasattr(bonded_atoms, '__iter__'):  # An iterable like a list
            for bonded_atom in bonded_atoms:
                atom.add_to_topology(bonded_atom)

    def del_atom(self, atom, delete_topology=True):
        """Deletes the specified :obj:`atom` from the molecule.

        Parameters
        ----------
        atom: :obj:`atom`
            The atom object to delete from the molcule
        delete_topology: bool
            If True, this atom will be removed from its bonded atom topologies.


        .. note:: This function invalidates cache objects with the attribute
                  or method 'preserve_cache_del_atoms'

        .. note:: Since atom objects are only referenced in residues, this
                  should delete the only atom reference in the molecule.
                  However, if the atom was referenced elsewhere, it won't get
                  deleted until that reference is deleted. An implementation
                  that returns weakref.proxy links from the molecule may be
                  preferable.

        .. note:: The delete_topology function will not remove hydrogens from
                  the atom's topology (though the H atom will not longer be
                  list in the atom's bonded_atoms). This is by design.

        Examples
        --------
        >>> mol = Molecule('2KXA')
        >>> print('{:.2f}, {:.2f}'.format(mol.mass, mol['A'][3].mass))
        2445.07, 147.19
        >>> sorted(mol['A'][3]['CA'].bonded_atoms())
        [A.F3.C, A.F3.CB, A.F3.HA, A.F3.N]
        >>> mol.del_atom(mol['A'][3]['N'])
        >>> print('{:.2f}, {:.2f}'.format(mol.mass, mol['A'][3].mass))
        2431.06, 133.18
        >>> sorted(mol['A'][3]['CA'].bonded_atoms())
        [A.F3.C, A.F3.CB, A.F3.HA]
        """
        # Invalidate the caches
        self.clear_cache(scope='del_atoms')

        if delete_topology:
            for bonded_atom in atom.bonded_atoms():
                bonded_atom.del_from_topology(atom)
        del atom.residue[atom.name]

    def strip_atoms(self, element):
        """Deletes all atoms with the given element (name, i.e. 'H')


        .. note:: This function invalidates cache objects with the attribute
                  or method 'preserve_cache_del_atoms'

        Examples
        --------
        >>> mol = Molecule('2KXA')
        >>> print('{:.2f}'.format(mol.mass))
        2445.07
        >>> mol.strip_atoms('H')
        >>> print('{:.2f}'.format(mol.mass))
        2285.49
        """
        # Invalidate the caches
        self.clear_cache(scope='del_atoms')

        # TODO: add option to also remove from topologies
        # TODO: Support '|' OR elements.
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
        first and last residues of each chain."""
  
        for chain in self.chains:
            # Get the residue numbers
            residues = list(chain.residues)
            no_residues = len(residues)

            # If it's a HETATM chain, every residue is it's own molecule,
            # and therefore they're each the first and last residue
            if '*' in chain.id:
                for residue in residues:
                    residue.first = True
                    residue.last = True
                    residue.prev_residue = None
                    residue.next_residue = None
                continue

            # Otherwise, create the linked list for each chain.
            # The first residue creates a new chain
            if no_residues > 0:
                first_residue = residues[0]
                first_residue.first = True
                first_residue.prev_residue = None

                if no_residues > 1:
                    second_residue = residues[1]
                    first_residue.last = False
                    first_residue.next_residue = second_residue
                else:
                    first_residue.last = True
                    first_residue.next_residue = None

            # Treat the last residue as special
            if no_residues > 0:
                last_residue = residues[-1]
                last_residue.last = True
                last_residue.next_residue = None

                if no_residues > 1:
                    before_last = residues[-2]
                    last_residue.first = False
                    last_residue.prev_residue = before_last
                else:
                    last_residue.first = True
                    last_residue.prev_residue = None

            # Now work on all the residues in between, starting with residue 2
            for count, residue in enumerate(residues[1:-1], 1):
                residue.first = False
                residue.last = False
                try:
                    residue.prev_residue = residues[count - 1]
                    residue.next_residue = residues[count + 1]
                except KeyError:
                    continue

    def set_atom_topologies(self):
        """Sets special topological information for specific atoms, like the
        amino terminal N or cystein bridges."""
        # TODO: Handle salt bridges
        # TODO: Handle hydrogen bonds
        # TODO: Handle XPLOR formatted H names like 'HN'

        # Properly set the topology for the alpha-amino and c-terminal groups
        # of each chain
        for chain in self.chains:
            first_residue = next(i for i in chain.residues if i.first)
            N = first_residue.get('N', None)
            if N is not None:
                if N.residue.name == 'PRO':
                    atom_list = ('H1', 'H2',)
                else:
                    atom_list = ('H1', 'H2', 'H3')

                for a in atom_list:
                    N.topology.update([a,])
                if 'HN' in N.topology:
                    N.topology.remove('HN')
                if 'H' in N.topology:
                    N.topology.remove('H')
                if 'C-1' in N.topology:
                    N.topology.remove('C-1')

            last_residue = next(i for i in chain.residues_reversed if i.last)
            C = last_residue.get('C', None)
            if C is not None:
                C.topology.update(['OXT',])
                if 'N+1' in C.topology:
                    C.topology.remove('N+1')
            OXT = last_residue.get('OXT', None)
            if OXT is not None:
                OXT.topology.update(['C','HXT'])

        # This part sets the connectivities listed under self.connections
        # Collect atom numbers and their respective atoms
        atom_numbers = []
        for connection in self.connections:
            # PDB CONECT fields have 11 fields
            if not len(connection) == 11:
                continue

            # 0: target_atom, 1: bonded_atom1, 2: bonded_atom2, 3: bonded_atom3
            # 4: bonded_atom4
            atom_numbers.append(connection[0:5])

        # Convert the atom numbers to actual atoms. First collect all the
        # atom numbers needed and convert it to a dict.
        number_set = set(ichain(*atom_numbers))
        atom_list = [a for a in self.atoms if a.number in number_set]
        atom_dict = {a.number:a for a in atom_list}

        # Set topological information
        for tgt_no, b1_no, b2_no, b3_no, b4_no in atom_numbers:
            tgt_atom = atom_dict.get(tgt_no, None)

            if tgt_atom is None:  # atom not found
                continue

            for b_no in (b1_no, b2_no, b3_no, b4_no):
                # Get the bonded atom's name a figure out if it's in the
                # same residue
                bonded_atom = atom_dict.get(b_no, None)
                if bonded_atom is None:  # atom not found
                    continue

                # add the atoms to each other topology.
                # Some groups are oxidized when forming a bond, like Cys
                # bridges, so a protons has to be removed from each atom in the
                # bond. This is accomplished by 'replace_in_topology;
                if tgt_atom != bonded_atom:
                    if (tgt_atom.residue.name, tgt_atom.name)  == ('CYS','SG'):
                        tgt_atom.replace_in_topology(bonded_atom, 'H')
                    else:
                        tgt_atom.add_to_topology(bonded_atom)

    def renumber_atoms(self):
        """Reset the atom numbers.

        This function is useful for resetting the atom numbers when atoms
        are removed or added.

        Parameters
        ----------
        skip_for_TER: bool (optional)

        Since the self.connections list depends on
        atom numbers, it is reset.


        .. note:: This function invalidates caches for objects without the
                  attribute or key 'preserve_cache_renumber_atoms'.
        """
        # Invalidate cache
        self.clear_cache(scope='renumber_atoms')

        self.connections = []

        # First order by the chains. The protein chains ('A', 'B', ...) should
        # come before the heteroatom chains ('A*', 'B*', ...)
        chains = list(self.chains)
        chains_protein = [c for c in sorted(chains, key = lambda c : c.id)
                          if not c.id.endswith('*')]
        chains_hetatm = [c for c in sorted(chains, key = lambda c : c.id)
                         if c.id.endswith('*')]
        chains = chains_protein + chains_hetatm

        # Prepare the atom number counter and renumber the atom numbers
        counter = count(1)
        for chain in chains:
            for residue in chain.residues:
                # The following sorts hydrogens together with their heavy atoms.
                # ex: ['N', 'H', 'C', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3']
                for atom in sorted(residue.atoms,
                                   key = lambda a : (a.name[1:], a.name[0])):
                    atom.number = next(counter)

    def clear_cache(self, scope=None):
        """Clear objects from this molecule's cache.

        Parameters
        ----------
        scope: str or None
            If None, all cache items are cleared.
            If a scope string is specified, all cache items that either do
            not have a preserve_cache_scope attribute equal to True or a
            key with 'preserve_cache_scope` will be removed. The scope can be
            a variety of values.

            - ``wb_rotation``: Cache items that are invalidated by a whole-body
              (whole molecule) rotation. Must have the
              preserve_cache_wb_rotation attribute equal to True or the
              'preserve_cache_wb_rotation` key.
            - ``wb_translation``: Cache items that are invalidated by whole-
              body translations.
            - ``add_atoms``: Cache items that are invalidated by adding atoms
            - ``del_atoms``: Cache items that are invalidated by deleting atoms
            - ``renumber_atoms``: Cache items that are invalidated by
              renumbering items.

        Returns
        -------
        int:
            The number of cache items deleted

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2KXA')
        >>> class Object: pass
        >>> obj = Object()
        >>> obj.preserve_cache_wb_rotation = True
        >>> mol.cache['obj'] = obj
        >>> mol.clear_cache('wb_rotation')
        0
        >>> 'obj' in mol.cache
        True
        >>> obj.preserve_cache_wb_rotation = False
        >>> mol.clear_cache('wb_rotation')
        1
        >>> 'obj' in mol.cache
        False
        """
        count = 0
        keys_to_delete = []

        # Find all of the cache kinds that should be deleted
        for k,v in self.cache.items():
            # Construct the name of the preserve_cache attribute or key
            preserve_name = ('_'.join(('preserve_cache', scope)) if scope
                             else None)

            # Skip item if it matches the preserve_name from the scope
            if preserve_name:
                if hasattr(v, '__iter__') and preserve_name in v:
                    continue
                if getattr(v, preserve_name, False):
                    continue
            keys_to_delete.append(k)

        # Setup debug messages
        msg = "Clearing cache scope '{}'"

        # Go through the keys_to_delete and delete these items from the cache
        for k in keys_to_delete:
            count += 1
            del self.cache[k]
            logging.debug(msg.format(scope if scope is not None else 'all'))

        return count

    # Read and Write Methods

    def write_pdb(self, filename):
        """Write data to a PDB file.

        Parameters
        ----------
        filename : str
            The path and filename to write the file to.
        """
        # FIXME: Add TER lines
        self.renumber_atoms()
        with open(filename, 'w') as f:
            # Populate 'ATOM' lines.
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
            # Populate 'CONECT' lines
            conect = set()
            conect_line = ("CONECT"
                           "{:>5}"
                           "{:>5}"
                           "{:>5}"
                           "{:>5}"
                           "\n")

            for atom in self.atoms:
                # If the atom is in the standard topology, then it is a
                # regular ATOM. Otherwise, treat it as a HETATM
                line_type = ('ATOM'
                             if atom.residue.name in topology else 'HETATM')

                # Prepare and write the 'ATOM' line
                atom_parms = {'line_type': line_type,
                              'atom_num': atom.number,
                              'atom_name': atom.name,
                              'alt_loc': '',
                              'res_name': atom.residue.name,
                              'chain': atom.chain.id[:1],  # remove * from ids
                              'res_number': atom.residue.number,
                              'icode': '',
                              'x': atom.pos[0],
                              'y': atom.pos[1],
                              'z': atom.pos[2],
                              'occupancy': 1,
                              'B_factor': 0,
                              'element': atom.element,
                              'charge': (atom.charge if hasattr(atom, 'charge')
                                         else '')}

                f.write(atom_line.format(**atom_parms))

                # Prepare the CONECT lines. CONECT lines aren't created for
                # Hydrogen atoms
                if atom.element == 'H' or atom.element == 'D':
                    continue

                if '*' in atom.chain.id:
                    # HETATM molecule chains are identified by a '*' in the
                    # chain ID. For these, all bonded atoms have to be
                    # saved in the CONECT lines.
                    bonded_atoms = atom.bonded_atoms()
                else:
                    # Other ATOMS whose topologies are known do not have a
                    # '*' in the chain ids. We just need to create CONECT
                    # lines for connectivities between residues
                    bonded_atoms = atom.bonded_heavy_atoms(longrange=True)
                bonded_numbers = sorted([b.number for b in bonded_atoms])

                if bonded_numbers:
                    for bonded_number in bonded_numbers:
                        # Skip duplicate connectivities
                        if (bonded_number, atom.number) in conect:
                            continue
                        conect.add((atom.number, bonded_number))

            # First the conect set has to be reformatted to group entries
            # for the same atom number.
            # FIXME: CONECT records shouldn't be made for Hydrogens, but heavy
            # atoms should have hydrogens in their CONECT records.
            conect_dict = dict()
            for atom_number_1, atom_number_2 in sorted(conect):
                atom_1_set = conect_dict.setdefault(atom_number_1, set())
                atom_1_set.add(atom_number_2)

            # Write the CONECT lines
            for atom_number_1, bonded_numbers in sorted(conect_dict.items()):
                bonded_numbers = sorted(bonded_numbers)

                # Each CONECT line only supports 3 bonded atoms, so a
                # separate CONECT LINE has to be added when there are more
                # than 3 bonded atoms. The grouper produces groups of 3 atom
                # numbers
                for numbers in grouper(3, bonded_numbers, ''):
                    f.write(conect_line.format(atom_number_1, *numbers))

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
        """Reads in data from a PDB file.

        Supported extensions include '.pdb' and '.pdb.gz'


        .. note:: This function invalidates all caches.
        """
        # Invalidate all caches
        self.clear_cache()

        if filename.endswith('.gz'):
            with gzip.open(filename) as f:
                self.read_stream(f)
        else:
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
        """
        temp_path = get_or_fetch(pdb_code, extensions='pdb.gz',
                                 urls=settings.pdb_urls)
        if temp_path is None:
            msg = "The identifier or file '{}' could not be found."
            raise IOError(msg.format(pdb_code))
        self.read_pdb(temp_path)

    _re_atom = re.compile((r"(?P<type>ATOM  |HETATM)"
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

    _conversions = {'type': lambda x: str(x).strip(),
                    'number': int,
                    'name': lambda x: str(x).strip(),
                    'alt_loc': lambda x: str(x).strip(),
                    'residue_name': lambda x: str(x).strip(),
                    'chain': lambda x: str(x).strip(),
                    'residue_number': int,
                    'icode': lambda x: str(x).strip(),
                    'x': float,
                    'y': float,
                    'z': float,
                    'occupancy': float,
                    'B_factor': float,
                    'element': lambda x: str(x).strip(),
                    'charge': lambda x: float(x) if str(x).strip() else ''}

    def _match_atom(self, match):
        """Matches an ATOM or HETATM line in a PDB file and creates new
        :obj:`mollib.Atom`, :obj:`mollib.Residue` and :obj:`mollib.Chain`
        objects.

        Parameters
        ----------
        match : regex :obj:`match`
            A regex match object generated by self._re_atom

        .. note:: This function uses the string conversion functions in
                  self._conversions.
        """
        t = match.groups()
        # Convert types
        t = (str(t[0]).strip(),  # 0: type
             int(t[1]),     # 1: number
             str(t[2]).strip(),  # 2: name
             str(t[3]).strip(),  # 3: alt_loc
             str(t[4]).strip(),  # 4: residue_name
             str(t[5]).strip(),  # 5: chain
             int(t[6]),     # 6: residue_number
             str(t[7]).strip(),  # 7: icode
             float(t[8]),   # 8: x
             float(t[9]),   # 9: y
             float(t[10]),  # 10: z
             float(t[11]),  # 11: occupancy
             float(t[12]),  # 12: B-factor
             str(t[13]).strip(), # 13: element
             str(t[14]).strip(), # 14: charge
             #(float(t[14]) if t[14] and
             # str(t[14]).strip() else '')  # 14: charge
             )

        # Implementation 2: twice as slow as implementation 1
        # groupdict = {field_name: self._conversions[field_name](field_value)
        #              for field_name, field_value
        #              in match.groupdict().items()}
        #
        # Implementation 3: slower than implementaiton 2
        # Old implementation of string conversion using the slower convert(...)
        # groupdict = {field_name: convert(field_value)
        #              for field_name, field_value
        #              in match.groupdict().items()}

        # create Chain, if it doesn't already exist
        identifier = t[5]

        # If this is a HETATM, then append a '*' to the chain name so that
        # it doesn't overwrite protein chains.
        if t[0] == 'HETATM':
            identifier += '*'

        # Create a new chain, if it doesn't already exist
        if identifier not in self:
            chain = self.chain_class(identifier=identifier)
            chain.molecule = self
            self[identifier] = chain

        chain = self[identifier]

        # create Residue, if it doesn't already exist
        res_number, res_name, = t[6], t[4]

        if res_number not in chain:
            try:
                residue = self.residue_class(number=res_number, name=res_name)

                residue.chain = chain
                residue.molecule = self
                chain[res_number] = residue
            except KeyError:
                return None
        residue = chain[res_number]

        # create the Atom. The following code overwrites atoms duplicate
        # in atom name
        number, name, alt_loc, element = t[1], t[2], t[3], t[13]

        # Reformat the x/y/z coordinates to a numpy array
        pos = np.array(t[8:11])

        # Populate the new atom parameters and create the new Atom object.
        atom = self.atom_class(number=number, name=name, pos=pos,
                               element=element, residue=residue, chain=chain,
                               molecule=self)

        residue[name] = atom

    _re_conect = re.compile(r'^CONECT(?P<numbers>[\s\d]+)$')

    def _match_conect(self, match):
        """Matches a 'CONECT' line in a PDB file and populates the connections
        attributue.
        """
        # The CONECT line format uses fixed columns and some may be missing so
        # This line has to be broken into fixed pieces
        number_str = match.groupdict()['numbers']

        str_len = len(number_str)
        offset = -7 # This is for the stripped 'CONECT' string
        cols = [(7, 11), (12, 16), (17, 21), (22, 26), (27, 31), (32, 36),
                (37, 41), (42, 46), (47, 51), (52, 56), (57,61)]
        cols = [(i+offset, j+offset+1) for i,j in cols]

        atom_numbers = [number_str[i:j].strip()
                        if (i < str_len and j <str_len) else None
                        for i,j in cols]
        atom_numbers = [int(i) if i else None for i in atom_numbers]
        self.connections.append(atom_numbers)

    def read_stream(self, stream):
        """Reads in data from a stream.
        """
        # TODO: Make a faster reader with structs and binary, with a fallback
        #       to the regex matcher.
        self.clear()

        # A list of regex matchers to harvest data from each line.
        # The first item is the regex to match. If there is a match,
        # the second item (function) will be called with the match
        matchers = OrderedDict((
                            ('ATOM', (self._re_atom, self._match_atom)),
                            ('CONECT', (self._re_conect, self._match_conect)),
                                ))

        # Find the ATOM/HETATM lines and pull out the necessary data

        # Generator implementation 1: This generator function reads only the
        # the first model. It's a little slower than implementation 2 (about
        # 10%). It takes 6.0s to read 3H0G
        for line in stream:
            m, func = None, None
            # Gzipped files return bytes lines that have to be decoded
            if type(line) == bytes:
                line = line.decode('latin-1')

            # Read only the first model. Skip reading ATOM lines hereafter
            if line[0:6] == 'ENDMDL' and 'ATOM' in matchers:
                del matchers['ATOM']

            # Advance the iterator if 1) no regex has been specified
            # yet, 2) the current regex is returning no matches and
            # it has before
            for name, (regex, func) in matchers.items():
                m = regex.match(line)
                if m:
                    break
            # If line cannot be matched, skip this line.
            if m is None or func is None:
                continue
            func(m)

        # Create links and set the atom topologies
        self.link_residues()

        # TODO: This test should skip CONECT items if there are more than
        #       99,999 atoms
        self.set_atom_topologies()

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
        # TODO: turn into method that defaults to excluding water molecule
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

    @property
    def pH(self):
        "Return the pH of the sample this molecule is in."
        if 'pH' in self._parameters:
            return self._parameters['pH']
        return settings.default_pH

    @pH.setter
    def pH(self, value):
        "Set the value of the pH for the sample of this molecule."
        self._parameters['pH'] = value