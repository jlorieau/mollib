import itertools


class Chain(dict):
    """A chain in a molecule.

    Parameters
    ----------
    id: str
        The chain's id. ex: 'A'
    molecule: `mollib.molecule`
        The molecule the chain belongs to.
    """

    def __init__(self, identifier, *args, **kwargs):
        self.id = identifier
        super(Chain, self).__init__(*args, **kwargs)

    def __repr__(self):
        return "{}".format(str(self.id))

    @property
    def residues(self):
        """An iterator over all residues in this chain,
        sorted by residue number."""
        return (residue for residue in
                sorted(self.values(), key=lambda a: a.number))

    @property
    def residues_reversed(self):
        """An iterator over all residues in this chain, sorted in *reversed*
        residue number."""
        return (residue for residue in
                sorted(self.values(), key=lambda a: a.number, reverse=True))

    @property
    def residue_size(self):
        """The number of residues in this chain."""
        return len(list(self.residues))

    @property
    def atoms(self):
        """An iterator over all atoms in this chain,
        sorted by atom number."""
        return (a for a in sorted(itertools.chain(*[r.values()
                                                    for r in self.values()]),
                                  key=lambda a: a.number))

    @property
    def atom_size(self):
        """The number of atoms in this chain."""
        return len(list(self.atoms))

    @property
    def mass(self):
        """The mass of all the atoms in this chain.

        >>> from mollib import Molecule
        >>> mol = Molecule('3C9J')
        >>> print("{:.2f}".format(mol['A'].mass)) # Chain A mass
        2507.61
        >>> print("{:.2f}".format(mol.mass)) # Total molecule mass
        10164.55
        """
        return sum(a.mass for a in self.atoms)
