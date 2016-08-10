from .primitives import Primitive


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
                 'bonded_atom_names')
    optional = ('charge', 'residue', 'chain', 'molecule',
                'bonded_atom_names')
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
