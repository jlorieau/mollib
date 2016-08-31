import re

from .primitives import Primitive
from .topology import *


# A regex to match the atom fullname
re_atom = re.compile(r'\s*(?P<molecule>[A-Z0-9]{1,4})\.'
                     r'(?P<chain_id>[A-Z]{1,3})\.'
                     r'(?P<residue_letter>[A-Z])'
                     r'(?P<residue_number>\d+)'
                     r'-'
                     r'(?P<atom_name>\w+)\s*')


def sorted_atom_list(atom_seq):
    """Sort the atoms in the given sequence into a list by stereochemical
    priority.

    Parameters
    ----------
    atom_seq : sequence
        A sequence type of (unsorted) :obj:`atom` objects.

    Returns
    -------
    list
        A sorted list of atoms.


    .. note:: The current implementation sorts by mass (reverse order), and
              only looks at one level of bonded_atoms to sort atoms with the
              same mass. A better implementation would seach bonded (and
              bonded of bonded) through a robust recursive function.

    Examples
    --------
    >>> from mollib import Molecule
    >>> mol = Molecule('2MJB')
    >>> I3 = mol['A'][3]
    >>> sorted_atom_list(I3['CA'].bonded_atoms())
    [I3-N, I3-C, I3-CB, I3-HA]
    >>> R42 = mol['A'][42]
    >>> sorted_atom_list(R42['CZ'].bonded_atoms())
    [R42-NE, R42-NH2, R42-NH1]
    """
    # First find the masses of all atoms. These will be used to sort
    # the atoms
    masses = {a: [a.mass, ] for a in atom_seq}
    mass_list = [a.mass for a in atom_seq]

    # Then find all duplicate masses. These are equivalent and have to
    # be sorted based on their bonded atoms.
    duplicate_masses = {k: v for k, v in masses.items()
                        if mass_list.count(v[0]) > 1}

    # No duplicates found. Just return the sorted masses
    if not duplicate_masses:
        return [a for a, m in sorted(masses.items(), reverse=True,
                                     key=lambda t: t[1])]

    # Otherwise look at the bonded atoms for the masses
    mass_of_bonded = {a: sum([b.mass for b in a.bonded_atoms()])
                      for a in duplicate_masses.keys()}

    # This loop extends the mass list in masses so that the bonded mass
    # can be used as a secondary sorting method. If these are still equal,
    # just use the atom's name
    for a, l in masses.items():
        if a in mass_of_bonded:
            l.append(mass_of_bonded[a])
            l.append(a.name)
        else:
            l.append(0)
            l.append(a.name)

    return [a for a, m in sorted(masses.items(), reverse=True,
                                 key=lambda t: t[1])]


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


    .. note:: __slots__ are used instead of a dicts because this can save a lot
              memory when many are instantiated and have faster attribute
              access.
    """

    # These are the required field. 'pos' (position)is the coordinate position
    # of the atom, as a numpy array
    __slots__ = ('number', 'name', 'pos', 'charge', 'element',
                 'residue', 'chain', 'molecule',
                 '_pK', '_topology')
    optional = ('charge', 'residue', 'chain', 'molecule',
                '_pK', '_topology')
    # bonded_atom_names ' ['1N', '2C-1'

    # Atom molecular weights. These must be labeled according the a str.title()
    # function. eg. ZN becomes Zn.
    atom_Mw = {'H': 1.01, 'D': 2.00, 'C': 12.01, '13C': 13.00, 'N': 14.01,
               '15N': 15.00, 'O': 16.00, 'Na': 22.99, 'Mg': 24.31, 'P': 30.97,
               'S': 32.07, 'Cl': 35.45, 'Zn': 65.38, 'Br': 79.904}

    def __repr__(self):
        return "{}-{}".format(self.residue, self.name) if self.residue else \
            "{}".format(self.name)

    def __lt__(self, other):
        # TODO: add chain comparisons
        return self.fullname.__lt__(other.fullname)

    def __le__(self, other):
        return self.fullname.__le__(other.fullname)

    def __eq__(self, other):
        return self.fullname.__eq__(other.fullname)

    def __ne__(self, other):
        return self.fullname.__ne__(other.fullname)

    def __gt__(self, other):
        return self.fullname.__gt__(other.fullname)

    def __ge__(self, other):
        return self.fullname.__ge__(other.fullname)

    def __hash__(self):
        return self.fullname.__hash__()

    @property
    def fullname(self):
        """The full name of this atom, including the molecule and chain.

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2KXA')
        >>> mol['A'][16]['CA'].fullname
        '2KXA.A.G16-CA'
        """
        molecule_name = (self.molecule.name
                         if self.molecule is not None else '')
        chain_id = (self.chain.id
                    if self.chain is not None else '')
        return '.'.join((molecule_name, chain_id, self.__repr__()))

    @property
    def mass(self):
        """The mass of this atom, depending on its element."""
        return self.atom_Mw[self.element.title()]

    # Atom coordinate getters and setters
    @property
    def x(self):
        """The x-coordinate of this atom."""
        return self.pos[0]

    @x.setter
    def x(self, value):
        self.pos[0] = value

    @property
    def y(self):
        """The y-coordinate of this atom."""
        return self.pos[1]

    @y.setter
    def y(self, value):
        self.pos[1] = value

    @property
    def z(self):
        """The z-coordinate of this atom."""
        return self.pos[2]

    @z.setter
    def z(self, value):
        self.pos[2] = value

    @property
    def topology(self):
        """The complete list of the topology of other atom *names* bonded
        this atom.

        Returns
        -------
        bonded_list : list
            A list of the atoms names (str) of all :obj:`atoms` that may be
            bonded to this :obj:`atom`.


        .. note:: :obj:`atoms` from the preceding residue are terminated with
                  '-1' and :obj:`atoms` from the proceeding residue are
                  terminated with '+1'.

        .. note:: The topology is a set, and can be modified using the standard
                  set operations

        .. note:: The topology list contains atoms that may or may not be
                  actually bonded to this atom--most notable hydrogen atoms.
                  Use the bonded_atoms to get a set of actual atoms bonded to
                  this atom.

        .. note:: The first letter of each topology string item must be the
                  atom's element. (i.e. 'CA')

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2PTN')
        >>> C22 = mol['A'][22]
        >>> sorted(C22['N'].topology)
        ['C-1', 'CA', 'H']
        >>> sorted(C22['SG'].topology)  # disulfide bridge
        ['2PTN.A.C157-SG', 'CB']
        """
        # TODO: Add Molecule functionality for cystein bridges
        # TODO: Add Molecule functionality to set first and last atom.
        if hasattr(self, '_topology'):
            return self._topology
        try:
            self._topology = topology[self.residue.name][self.name].copy()
        except KeyError:
            self._topology = set()
        return self._topology

    @topology.setter
    def topology(self, value):
        self._topology = value

    def reset_topology(self):
        "Reset this atom's topology to its default values."
        if hasattr(self, '_topology'):
            del self._topology

    def add_to_topology(self, atom):
        """Adds the given atom to this atom's topology.


        .. note:: Since the topology is a set, there is no risk in adding
                  an atom to the topology twice.

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2PTN') # structure with no Hs
        >>> G18 = mol['A'][18]
        >>> sorted(G18['N'].topology)
        ['C-1', 'CA', 'H']
        >>> G18['N'].add_to_topology(G18['C'])
        >>> sorted(G18['N'].topology)
        ['C', 'C-1', 'CA', 'H']
        >>> G18['N'].add_to_topology(mol['A'][16]['N'])
        >>> sorted(G18['N'].topology)
        ['2PTN.A.I16-N', 'C', 'C-1', 'CA', 'H']
        """
        # If the residues are the same, simply add the atom's name
        if (self.residue is not None and atom.residue is not None and
            self.residue == atom.residue):
            self_name = self.name
            atom_name = atom.name
        # Otherwise use the atom's fullname
        else:
            self_name = self.fullname
            atom_name = atom.fullname

        # Topologies are sets, so there is no risk of adding duplication
        # atom names.
        self.topology.update([atom_name,])
        atom.topology.update([self_name,])

    def del_from_topology(self, atom):
        """Delete atom from this atom's topology and vice-versa.

        Parameters
        ----------
        atom: :obj:`atom`
            The atom to remove from this atom's topology.

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2KXA') # structure with no Hs
        >>> F3 = mol['A'][3]
        >>> sorted(F3['N'].topology)
        ['C-1', 'CA', 'H']
        >>> F3['N'].del_from_topology(F3['H'])
        >>> sorted(F3['N'].topology)
        ['C-1', 'CA']
        """
        # If the residues are the same, simply add the atom's name
        if (self.residue is not None and atom.residue is not None and
                    self.residue == atom.residue):
            self_name = self.name
            atom_name = atom.name
        # Otherwise use the atom's fullname
        else:
            self_name = self.fullname
            atom_name = atom.fullname

        self.topology -= {atom_name,}
        atom.topology -= {self_name,}

    def in_topology(self, atom):
        """Test whether atom is already in this atom's topology."""
        if (self.residue is not None and atom.residue is not None and
            self.residue.name == atom.residue.name):
            atom_name = atom.name

            # Do not set the _topology cache if not needed
            if not hasattr(self, '_topology'):
                try:
                    return atom_name in topology[self.residue.name][self.name]
                except KeyError:
                    return atom_name in self.topology

        # Otherwise use the atom's fullname
        else:
            atom_name = atom.fullname
            return atom_name in self.topology

    def replace_in_topology(self, atom, start_str='H'):
        """Replaces an atom name from the topology starting with start_str
        with the :obj:`atom` object.

        This function is useful for creating new bonds while removing default
        Hs bound to an atom, like a cystein SG.

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2KXA')
        >>> A5 = mol['A'][5]
        >>> sorted(A5['N'].topology)
        ['C-1', 'CA', 'H']
        >>> A5['N'].replace_in_topology(A5['CB'])
        >>> sorted(A5['N'].topology)
        ['C-1', 'CA', 'CB']
        >>> sorted(A5['CB'].topology) # both atoms lose an H.
        ['CA', 'HB2', 'HB3', 'N']
        """
        # Go through the topologies of each atom.
        for a1, a2 in ((self, atom), (atom, self)):
            # Find a an atom name to replace. Use the first found.
            try:
                # sort the topology set so that the proton remove is
                # deterministic
                remove = next(i for i in sorted(a1.topology)
                              if i.startswith(start_str))

                # Note that an atom should only be replaced if the the given
                # atom is not already in the topology.
                if not a1.in_topology(a2):
                    a1.topology.remove(remove)
            except StopIteration:
                continue

        self.add_to_topology(atom)

    def bonded_atoms(self, sorted=False):
        """The atoms bonded to this atom, based on the topology method.

        Parameters
        ----------
        sorted: bool
            If True, atoms will be sorted according to their stereochemical
            priority

        Returns
        -------
        bonded : set
          A set of the *actual* :obj:`atom` objects currently bonded to this
          atom.

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2PTN')
        >>> C22 = mol['A'][22]
        >>> C22['C'].bonded_atoms(sorted=True)
        [C22-O, G23-N, C22-CA]
        >>> C22['SG'].bonded_atoms(sorted=True)  # disulfide bridge
        [C157-SG, C22-CB]
        """
        bonded = set()
        for name in self.topology:
            # Retrieve an atom from the previous or next residue, if specified
            if name in self.residue:
                bonded |= {self.residue[name]}
            elif name.endswith('-1'):
                bonded |= {self.residue.prev_residue[name[:-2]]}
            elif name.endswith('+1'):
                bonded |= {self.residue.next_residue[name[:-2]]}
            # Retrieve an atom if it's using its fullname.
            # ex: 2PTN.A.C220-SG
            elif '-' in name:
                match = re_atom.match(name)
                if match:

                    residue_number = int(match.groupdict()['residue_number'])
                    atom_name = match.groupdict()['atom_name']
                    chain_id = match.groupdict()['chain_id']
                    molecule_name = match.groupdict()['molecule']

                    # FIXME: Current implementation only works for atoms in the
                    # same chain and molecule.
                    if (molecule_name != self.molecule.name or
                        chain_id != self.chain.id):
                        continue

                    residue = self.chain.get(residue_number, None)
                    atom = (residue.get(atom_name, None)
                            if residue is not None else None)
                    if atom is not None:
                        bonded |= {atom}
            else:
                continue
        if sorted:
            return sorted_atom_list(bonded)
        else:
            return list(bonded)

    def bonded_heavy_atoms(self, sorted=False):
        """The heavy atoms bonded to this atom, based on the topology method.

        Parameters
        ----------
        sorted: bool
            If True, atoms will be sorted according to their stereochemical
            priority

        Returns
        -------
        bonded : set
            A set of the *actual* :obj:`atom` objects for heavy atoms currently
            bonded to this atom.

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2PTN')
        >>> C22 = mol['A'][22]
        >>> C22['CA'].bonded_heavy_atoms(sorted=True)
        [C22-N, C22-CB, C22-C]
        """
        bonded = {a for a in self.bonded_atoms()
                        if not a.element == 'H' or a.element == 'D'}
        if sorted:
            return sorted_atom_list(bonded)
        else:
            return list(bonded)

    def bonded_interresidue_atoms(self, sorted=False):
        """The heavy atoms bonded on this atom that are from other residues.

        Parameters
        ----------
        sorted: bool
            If True, atoms will be sorted according to their stereochemical
            priority

        Returns
        -------
        bonded : set
            A set of the *actual* :obj:`atom` objects for atoms bonded to this
            atom from other residues.

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2PTN')
        >>> C22 = mol['A'][22]
        >>> C22['SG'].bonded_interresidue_atoms(sorted=True)
        [C157-SG]
        """
        bonded = set()

        # Atoms from other residues are only specified if this atom has
        # marked them in its topology. If the topology isn't set (ie:
        # _topology isn't set), then there are no atoms from other residues
        # bonded to this atom
        if not hasattr(self, '_topology'):
            return list(bonded)
        topology = self._topology

        # Atoms from other residues are identified by their fullname, which
        # include '.' characters
        for atom_name in topology:
            if '-' not in atom_name:
                continue

            match = re_atom.match(atom_name)
            if not match:
                continue

            molecule = match.groupdict()['molecule']
            if not self.molecule.name == molecule:
                continue
            molecule = self.molecule

            chain_id = match.groupdict()['chain_id']
            residue_number = int(match.groupdict()['residue_number'])
            name = match.groupdict()['atom_name']

            try:
                bonded_atom = molecule[chain_id][residue_number][name]
            except KeyError:
                continue

            bonded.add(bonded_atom)

        if sorted:
            return sorted_atom_list(bonded)
        else:
            return list(bonded)

    # @property
    # def geometry_atoms(self):
    #     """The other non-bonded atoms that impact the geometry of bonding at
    #     this atom.
    #
    #     For example, to determine the protonation of a C-O in a carboxylate,
    #     the position of the other 'O' will be needed.
    #
    #     Returns
    #     -------
    #     geometry_list : list
    #         A list of the other, *non-bonded* :obj:`atom` objects that
    #         influence this :obj:`atom`'s geometry.
    #
    #     Examples
    #     --------
    #     >>> from mollib import Molecule
    #     >>> mol = Molecule('2KXA')
    #     >>> D19 = mol['A'][19]
    #     >>> D19['OD1'].bonded_atoms
    #     ['D19-OD2']
    #     """
    #     raise NotImplementedError
    #
    # @property
    # def configuration(self):
    #     """
    #
    #     Returns
    #     -------
    #
    #     Examples
    #     --------
    #     >>> from mollib import Molecule
    #     >>> mol = Molecule('2KXA')
    #     >>> G1 = mol['A'][1]
    #     >>> mol.pH = 2.0
    #     >>> G1['N'].configuration
    #     'tetrahedral'
    #     >>> mol.pH = 12.0
    #     >>> G1['N'].configuration
    #     'tetragonal'
    #
    #
    #     .. note:: this only works for atoms without an ionization_alternative.
    #               It requires the pK attribute set and the molecule.pH
    #               attribute set.
    #     """
    #     raise NotImplementedError

