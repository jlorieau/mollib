import re

from .primitives import Primitive
from .topology import *


#: The regex to match the atom fullname.
re_atom_str = (r'\s*((?P<molecule>[A-Z0-9]{1,4})\:)?'  # Optional
               r'((?P<chain_id>[A-Z]{1,3})\.)?'        # Optional
               r'(?P<residue_letter>[A-Z])?'           # Optional
               r'(?P<residue_number>\d+)'
               r'\.?'                                  # Optional
               r'(?P<atom_name>[A-Z]+\d*)\s*')
# re_atom_str = (r'^\s*((?P<molecule>[A-Z0-9]{1,4})\:)?'  # Optional
#                r'((?P<chain_id>[A-Z]{1,3})\.)?'        # Optional
#                r'(?P<residue_letter>[A-Z])?'           # Optional
#                r'(?P<residue_number>\d+)'
#                r'\.?'                                  # Optional
#                r'(?P<atom_name>\w+)\s*$')

re_atom = re.compile(re_atom_str)


def sorted_atom_list(atom_seq):
    """Sort the atoms in the given sequence into a list by stereochemical
    priority.

    Parameters
    ----------
    atom_seq : sequence
        A sequence type of (unsorted) :obj:`Atom` objects.

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
    [A.I3.N, A.I3.C, A.I3.CB, A.I3.HA]
    >>> R42 = mol['A'][42]
    >>> sorted_atom_list(R42['CZ'].bonded_atoms())
    [A.R42.NE, A.R42.NH2, A.R42.NH1]
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

    # Otherwise look at the bonded atoms for the masses. Get the sum of these
    # masses. The sum of masses have to be rounded to avoid non-deterministic
    # floating point sums that change the order, like: 25.04 vs 25.040000001
    # one time and 25.0400000001 vs 25.04 the next.
    mass_of_bonded = {a: round(sum([b.mass for b in a.bonded_atoms()]), 2)
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

    .. note:: Atom objects support array access to the pos `numpy.array`.
              Items 0, 1, 2 correspond to the pos[0], pos[1] and pos[2].
              
    .. note:: Atom labels can take on a variety of values. The molecule name,
              chain id and residue letter are optional. The residue number and
              atom name are required.
              
              The following are valid labels:
                  - 2KXA:A.G16.HA2
                  - 2KXA:A.G16HA2
                  - A.G16HA2
                  - 16HA2
    """

    # These are the required field. 'pos' (position)is the coordinate position
    # of the atom, as a numpy array
    __slots__ = ('number', 'name', 'pos', 'charge', 'element',
                 'residue', 'chain', 'molecule',
                 '_pK', '_topology','_repr', '_fullname')
    optional = ('charge', 'residue', 'chain', 'molecule',
                '_pK', '_topology', '_repr', '_fullname')

    # Atom molecular weights. These must be labeled according the a str.title()
    # function. eg. ZN becomes Zn.
    atom_Mw = {'H': 1.01, 'D': 2.00, 'Be': 9.01, 'B': 10.81, 'C': 12.01,
               '13C': 13.00, 'N': 14.01, '15N': 15.00, 'O': 16.00, 'F':19.00,
               'Na': 22.99, 'Mg': 24.31, 'P': 30.97, 'S': 32.07, 'Cl': 35.45,
               'V': 50.94, 'Cr': 52.00, 'Fe': 55.84,'Ni': 58.69, 'Co': 58.93,
               'Cu': 63.55, 'Zn': 65.38, 'Ga': 69.72, 'Se': 78.96, 'As': 74.92,
               'Br': 79.904, 'Zr': 91.22, 'Mo': 95.94, 'Ru': 101.07,
               'I': 126.90, 'Te': 127.60, 'W': 183.84, 'Ir': 192.22,
               'Pt': 195.08, 'Hg': 200.59, 'Pb': 207.2, 'U': 238.05}

    def __repr__(self):
        if hasattr(self, '_repr'):
            return self._repr
        repr = ''
        if self.chain:
            repr += "{}.".format(self.chain)
        if self.residue:
            repr += "{}.".format(self.residue)
        repr += self.name
        self._repr = repr
        return repr

    def __lt__(self, other):
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

    def __getitem__(self, b):
        return self.pos[b]

    def __setitem__(self, b, value):
        self.pos[b] = value

    def __len__(self):
        return len(self.pos)

    @property
    def fullname(self):
        """The full name of this atom, including the molecule and chain.

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2KXA')
        >>> mol['A'][16]['CA'].fullname
        '2KXA:A.G16.CA'
        """
        if not hasattr(self, '_fullname'):
            molecule_name = (self.molecule.name
                             if self.molecule is not None else '')
            self._fullname = ':'.join((molecule_name, self.__repr__()))
        return self._fullname

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


        .. note:: :obj:`Atom` objects from the preceding residue are terminated 
                  with '-1' and :obj:`atoms` from the proceeding residue are
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
        ['2PTN:A.C157.SG', 'CB']
        """
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
        ['2PTN:A.I16.N', 'C', 'C-1', 'CA', 'H']
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
        """Test whether atom is already in this atom's topology.

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2KXA')
        >>> L2, F3, G4 = mol['A'][2], mol['A'][3], mol['A'][4]
        >>> F3['N'].in_topology(F3['CA'])
        True
        >>> F3['N'].in_topology(L2['C'])
        True
        >>> F3['C'].in_topology(G4['N'])
        True
        >>> F3['C'].in_topology(F3['N'])
        False
        """
        # Get the topology dict. Do not set the _topology cache if not needed
        top = (self._topology
               if hasattr(self, '_topology')
               else topology.get(self.residue.name, {}).get(self.name, ''))

        if self.chain.id == atom.chain.id:
            # If the two atoms are in residues that are separated by +/- 1 or 0,
            # then the atom name is just the Atom's name with +1, -1 or nothing
            # on the end.
            delta = atom.residue.number - self.residue.number
            atom_name = ("{name}{delta:+}".format(name=atom.name, delta=delta)
                         if delta != 0 else atom.name)

            # See if the atom's name is in there
            if atom_name in top:
                return True

        # Otherwise try the atom's full name. This can happen for covalent
        # bonds across the molecule, like Cystein bridges
        atom_name = atom.fullname
        return atom_name in top

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

    def bonded_atoms(self, sorted=False, longrange=False):
        """The atoms bonded to this atom, based on the topology method.

        Parameters
        ----------
        sorted: bool (optional)
            If True, atoms will be sorted according to their stereochemical
            priority
        longrange: bool (optional)
            If True, only atoms from other residues will be returned (not
            including residues +/- 1)

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
        [A.C22.O, A.G23.N, A.C22.CA]
        >>> C22['SG'].bonded_atoms(sorted=True)  # disulfide bridge
        [A.C157.SG, A.C22.CB]
        >>> C22['SG'].bonded_atoms(sorted=True, longrange=True)
        [A.C157.SG]
        """
        bonded = set()
        for name in self.topology:
            # Retrieve an atom from the previous or next residue, if specified.
            # But only do this if longrange is False
            if name in self.residue and not longrange:
                bonded |= {self.residue[name]}
                continue
            elif name.endswith('-1') and not longrange:
                try:
                    bonded |= {self.residue.prev_residue[name[:-2]]}
                except (KeyError, TypeError):
                    pass
                continue
            elif name.endswith('+1') and not longrange:
                try:
                    bonded |= {self.residue.next_residue[name[:-2]]}
                except (KeyError, TypeError):
                    pass
                continue

            # Retrieve an atom if it's using its fullname.
            # ex: 2PTN.A.C220-SG
            match = re_atom.match(name)
            if match:
                residue_number = int(match.groupdict()['residue_number'])
                atom_name = match.groupdict()['atom_name']
                chain_id = match.groupdict()['chain_id']
                molecule_name = match.groupdict()['molecule']

                # FIXME: Current implementation only works for atoms in the
                # same  molecule.
                if molecule_name != self.molecule.name:
                    continue

                chain = self.molecule.get(chain_id, None)

                residue = (chain.get(residue_number, None)
                           if chain is not None else None)

                atom = (residue.get(atom_name, None)
                        if residue is not None else None)

                # Skip non-longrange atoms, if specified
                if longrange and (self.chain != chain or
                                  abs(self.residue.number -
                                      residue_number) <= 1):
                    continue

                if atom is not None:
                    bonded |= {atom}
            else:
                continue
        if sorted:
            return sorted_atom_list(bonded)
        else:
            return list(bonded)

    def bonded_heavy_atoms(self, sorted=False, longrange=False):
        """The heavy atoms bonded to this atom, based on the topology method.

        Parameters
        ----------
        sorted: bool (optional)
            If True, atoms will be sorted according to their stereochemical
            priority
        longrange: bool (optional)
            If True, only atoms from other residues will be returned (not
            including residues +/- 1)

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
        [A.C22.N, A.C22.CB, A.C22.C]
        >>> C22['SG'].bonded_heavy_atoms(longrange=True)
        [A.C157.SG]
        """
        bonded = [a for a in self.bonded_atoms(sorted, longrange)
                  if not a.element == 'H' or a.element == 'D']
        return list(bonded)
