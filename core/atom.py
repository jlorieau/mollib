from .primitives import Primitive
from .topology import *
from . import settings


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
    atom_Mw = {'H': 1.01, 'C': 12.01, '13C': 13.00, 'N': 14.01, '15N': 15.00,
               'O': 16.00, 'Na': 22.99, 'Mg': 24.31, 'P': 30.97, 'S': 32.07,
               'Cl': 35.45, 'Zn': 65.38, }

    def __repr__(self):
        return u"{}-{}".format(self.residue, self.name) if self.residue else \
            u"{}".format(self.name)

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
        """The complete list of the topology of other atom *names* to bonded
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

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2KXA')
        >>> D19 = mol['A'][19]
        >>> D19['OD1'].topology
        {'CG'}
        >>> D19['OD1'].topology |= {'HD1'}
        >>> D19['OD1'].topology
        {'CG', 'HD1'}
        """
        if hasattr(self, '_topology'):
            return getattr(self, '_topology')
        try:
            return topology[self.residue.name][self.name]
        except KeyError:
            return set()

    @topology.setter
    def topology(self, value):
        if not hasattr(self, '_topology'):
            self._topology = self.topology
        self._topology = value

    # @property
    # def bonded_atoms(self):
    #     """The atoms bonded to this atom.
    #
    #     Returns
    #     -------
    #     bonded : list
    #       A list of the *actual* :obj:`atom` objects currently bonded to this
    #       atom.
    #
    #
    #     .. note:: This list may not correspond to the bonded_topology when
    #               there are missing :obj:`atom` objects, like hydrogens.
    #
    #     Examples
    #     --------
    #     >>> from mollib import Molecule
    #     >>> mol = Molecule('2KXA')
    #     >>> D19 = mol['A'][19]
    #     >>> D19['C'].bonded_atoms
    #     ['D19-CA', 'G20-N']
    #     """
    #     raise NotImplementedError
    #
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

