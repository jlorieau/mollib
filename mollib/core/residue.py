from collections import namedtuple

from .geometry import measure_dihedral
from . import settings

class IonGroup(namedtuple('IonGroup', ['possible_atoms', 'pKs'])):
    """A group of ionizeable atoms and their pKs.

    Attributes
    ----------
    possible_atoms: list of :obj:`mollib.Atom` objects
        A list of atoms that may be ionized. (i.e. may receive or lose a
        hydrogen atom, depending on pH)
    pKs: tuple of float
        A tuple with all of the pK values for the :obj:`mollib.Atom` objects
        in possible atoms.
    """
    pass


class ResidueProcessingError(Exception):
    "An error in processing the residue."
    pass


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
    chain: :obj:`mollib.chain`
        The chain the residue belongs to.
    molecule: :obj:`mollib.molecule`
        The molecule the residue belongs to.

    Attributes
    ----------
    one_letter_codes: dict
        A dict for converting 3-letter residue names (key) to 1-letter residue
        names (value). The residue object will use 'X' if the 3-letter name
        is not found.
    prev_residue: :obj:`Residue`
        The preceding residue object, in the sequence, from the :obj:`chain`.
        This is a link in a doubly-linked list.
    next_residue: :obj:`Residue`
        The proceeding residue object, in the sequence, from the :obj:`chain`.
        This is a link in a doubly-linked list.
    first: bool
        True if this is the first :obj:`residue` in the :obj:`chain`.
    last: bool
        True if this is the last :obj:`residue` in the :obj:`chain`.

        .. note:: The `prev_residue`, `next_residue`, `first` and `last` are
                  populated by the :obj:`molecule` object on creation--not the
                  :obj:`residue` object.
    """

    # These are set by the molecule
    molecule = None
    chain = None

    first = False
    last = False

    one_letter_codes = {'ALA': 'A', 'GLY': 'G', 'SER': 'S', 'THR': 'T',
                        'MET': 'M', 'CYS': 'C', 'ILE': 'I', 'LEU': 'L',
                        'VAL': 'V', 'PHE': 'F', 'TYR': 'Y', 'TRP': 'W',
                        'ASN': 'N', 'GLN': 'Q', 'ASP': 'D', 'GLU': 'E',
                        'HIS': 'H', 'PRO': 'P', 'ARG': 'R', 'LYS': 'K'}

    def __init__(self, name, number, *args, **kwargs):
        name = str(name).upper()

        self.molecule = None
        self.chain = None

        self.name = name  # full name, MET
        self.letter = self.one_letter_codes.get(self.name, 'X')  # letter code
        self.number = number

        # Lazily calculated variables
        self._classification = None

        super(Residue, self).__init__(*args, **kwargs)

    def __repr__(self):
        return "{}{}".format(self.letter, self.number)

    @property
    def atoms(self):
        """An iterator over all atoms in this residue, sorted by atom number"""
        return (atom for atom in sorted(self.values(), key=lambda a: a.number))

    @property
    def atom_size(self):
        """The number of atoms in this residue."""
        return len(list(self.atoms))

    @property
    def mass(self):
        """The mass of all atoms in this residue.

        >>> from mollib import Molecule
        >>> mol = Molecule('2KXA')
        >>> print("{:.2f}".format(mol['A'][3].mass)) # Phe-3 mass
        147.19
        """
        return sum(a.mass for a in self.atoms)

    @property
    def next_residue(self):
        "Return the next residue in the molecule, or None if not available."
        if hasattr(self, '_next_residue_number'):
            return self.chain[self._next_residue_number]
        else:
            return None

    @next_residue.setter
    def next_residue(self, residue):
        "Set the next residue in the molecule."
        if residue is None:
            return None
        if residue.molecule != self.molecule:
            msg = 'Next residue {} has to be from the same molecule {}.'
            raise ResidueProcessingError(msg.format(residue, self.molecule))
        if residue.chain != self.chain:
            msg = 'Next residue {} has to be from the same chain {}.'
            raise ResidueProcessingError(msg.format(residue, self.chain))
        self._next_residue_number = residue.number

    @property
    def prev_residue(self):
        "Return the previous residue in the molecule, or None if not available."
        if (hasattr(self, '_prev_residue_number') and
                    self._prev_residue_number in self.chain):
            return self.chain[self._prev_residue_number]
        else:
            return None

    @prev_residue.setter
    def prev_residue(self, residue):
        "Set the previous residue in the molecule."
        if residue is None:
            return None
        if residue.molecule != self.molecule:
            msg = 'Next residue {} has to be from the same molecule {}.'
            raise ResidueProcessingError(msg.format(residue, self.molecule))
        if residue.chain != self.chain:
            msg = 'Next residue {} has to be from the same chain {}.'
            raise ResidueProcessingError(msg.format(residue, self.chain))
        self._prev_residue_number = residue.number

    @property
    def ramachandran_angles(self):
        """The backbone Ramachandran angles for this residue.

        Returns
        -------
        angles: list
            A pair of angles (floats) representing the phi and psi angles in
            degrees.

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2KXA')
        >>> print("{:.1f} {:.1f}".format(*mol['A'][3].ramachandran_angles))
        -65.2 -48.3
        >>> print("{} {:.1f}".format(*mol['A'][1].ramachandran_angles))
        None -179.9
        """
        # Retrieve the appropriate atoms to calculate the ramachandran angles
        c_prev = (self.prev_residue.get('C', None)
                  if self.prev_residue is not None else None)
        n = self.get('N', None)
        ca = self.get('CA', None)
        c = self.get('C', None)
        n_next = (self.next_residue.get('N', None)
                  if self.next_residue is not None else None)

        angles = [None, None]  # Initial phi/psi angles are set to None
        # Iterate over phi (angles[0]) and phi (angles[1]) to calculate the
        # dihedrals
        for item, a, b, c, d in [(0, c_prev, n, ca, c),   # phi
                                 (1, n, ca, c, n_next)]:  # psi

            if a is None or b is None or c is None or d is None:
                continue
            angles[item] = measure_dihedral(a, b, c, d)
        return angles

    @property
    def classification(self):
        """The secondary structure classification for this residue.
        
        The classification is a 2-member tuple of strings with the major
        classification (ex: 'alpha-helix') and the minor classification
        ('N-term')
        
        Returns
        -------
        tuple or None
            The classification tuple or,
            None, if a classification has not yet been assigned.
        """
        return self._classification

    @classification.setter
    def classification(self, value):
        assert isinstance(value, tuple) and len(value) > 0
        self._classification = value

    @property
    def ionizeable_groups(self):
        """The number and identity of ionizeable groups (:obj:`IonGroup`) for
         this residue.

         The returned list is *independent* of the molecule's pH. It's just a
         listing of atom's that can be ionized in this residue.

        Returns
        -------
        list of :obj:`IonGroup`
            The list of ionizeable groups (:obj:`IonGroup`) for this residue.

        Examples
        --------
        >>> from mollib import Molecule
        >>> mol = Molecule('2PTN')
        >>> I16 = mol['A'][16] # First residue
        >>> I16.ionizeable_groups
        [IonGroup(possible_atoms=[A.I16.N], pKs=(7.7, 14.0, 14.0))]
        >>> H40 = mol['A'][40]
        >>> H40.ionizeable_groups
        [IonGroup(possible_atoms=[A.H40.ND1, A.H40.NE2], pKs=(6.6, 14.0))]
        >>> N245 = mol['A'][245] # Last residue
        >>> N245.ionizeable_groups
        [IonGroup(possible_atoms=[A.N245.O, A.N245.OXT], pKs=(-1.0, 3.3))]
        """
        groups = []
        # Create an ionization group for the alpha-amino or terminal-carboxylate
        if self.first:
            groups += settings.pKs['first'].items()
        if self.last:
            groups += settings.pKs['last'].items()
        if self.name in settings.pKs:
            groups += settings.pKs[self.name].items()

        # Convert the pKs and names to actual atom objects
        return_list = []
        for names, pKs in groups:
            # Note that multiple atom names are split with a '-' character
            try:
                atoms = [self[name] for name in names.split('-')]
            # Skip this group if an atom isn't found
            except KeyError:
                continue
            iongroup = IonGroup(possible_atoms=atoms, pKs=pKs)
            return_list.append(iongroup)
        return return_list
