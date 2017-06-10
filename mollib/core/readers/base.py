"""
Molecular reader classes.
"""
from mollib.core import Molecule, Chain, Residue, Atom


class MoleculeReader(object):
    """Molecule Reader Factory base class.
    
    This class is the interface used by read molecules. Sub-classes implement
    the different readers.
    """

    #: The order to run the reader (int). To disable, set the order to None.
    order = None

    #: The default class to create molecules
    molecule_class = Molecule

    #: The default class to create chains
    chain_class = Chain

    #: The default class to create residues
    residue_class = Residue

    #: The default class to create atoms
    atom_class = Atom

    def __init__(self, identifiers_or_files, *args, **kwargs):
        """Initialize the molecular reader
        
        Parameters
        ----------
        identifiers_or_files: str or list
            One or more identifiers (ex: PDB codes, like '2kxa') or
            filenames and paths for structure files.
        args: tuple, optional
            Positional arguments
        kwargs: dict, optional
            Keywork arguments
        """
        self._identifiers_or_files = identifiers_or_files
        super(MoleculeReader, self).__init__(*args, **kwargs)

    def parse(self, stream, molecule=None):
        """Parse a data stream.
        
        Parameters
        ----------
        stream: the data stream to parse.
        
        Returns
        -------
        count: int
            The number of matched items.
        """
        return 0
