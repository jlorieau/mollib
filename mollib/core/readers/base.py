"""
Molecular reader classes.
"""
import os
import gzip
import io

from .. import molecule
from .. import chain
from .. import residue
from .. import atom
from mollib.utils.net import get_or_fetch
from mollib.utils.checks import check_file


class MoleculeReader(object):
    """Molecule Reader Factory base class.
    
    This class is the interface used by read molecules. Sub-classes implement
    the different readers.
    
    Examples
    --------
    >>> from mollib import MoleculeReader
    >>> mr = MoleculeReader()
    >>> molecules = mr.read('2kxa')  # Read in the 10 models
    >>> print(len(molecules))
    10
    >>> print(molecules[0])
    Molecule (2kxa-1):    1 chains, 24 residues, 332 atoms.
    >>> molecules = mr.read('2kxa', model_ids=[1, 5, 8])
    >>> for i in molecules:
    ...     print(i)
    Molecule (2kxa-1):    1 chains, 24 residues, 332 atoms.
    Molecule (2kxa-5):    1 chains, 24 residues, 332 atoms.
    Molecule (2kxa-8):    1 chains, 24 residues, 332 atoms.
    """

    #: The order to run the reader (int). To disable, set the order to None.
    order = None

    #: A tuple of urls to check for the file or identifier
    urls = None

    #: The default class to create molecules
    molecule_class = molecule.Molecule

    #: The default class to create chains
    chain_class = chain.Chain

    #: The default class to create residues
    residue_class = residue.Residue

    #: The default class to create atoms
    atom_class = atom.Atom

    def __init__(self, urls=None, *args, **kwargs):
        """Initialize the molecular reader
        
        Parameters
        ----------
        urls: tuple
            A list of urls to check for the file to read
        args: tuple, optional
            Positional arguments
        kwargs: dict, optional
            Keywork arguments
        """
        if urls is not None:
            self.urls = urls

        # Find subclasses and instantiate these.
        self.subclasses = [cls(urls, *args, **kwargs)
                           for cls in sorted(self.__class__.__subclasses__(),
                                             key=lambda x: x.order)
                           if cls.order is not None]

    def read(self, identifiers_or_files, source_molecules=None,
             *args, **kwargs):
        """Read data into molecules.
        
        Parameters
        ----------
        identifiers_or_files: str or list
            One or more identifiers (ex: PDB codes, like '2kxa') or
            filenames and paths for structure files.
        source_molecule: molecule (:obj:`mollib.Molecule`), optional
            If specified, the molecule will have their data replaced with the 
            parsed data. Only the first model will be loaded.
        args: tuple, optional
            Positional arguments
        kwargs: dict, optional
            Keywork arguments
            
        Returns
        -------
        molecules: list of :obj:`mollib.Molecule` objects
            The molecules read from the file.
        """
        # Wrap the identifers or files into lists, if needed
        if (not hasattr(identifiers_or_files, '__iter__')
            or isinstance(identifiers_or_files, str)):
            identifiers_or_files = [identifiers_or_files, ]

        # Prepare the molecules to return
        molecules = []

        # Loop through the identifiers and files
        for identifier in identifiers_or_files:
            # The following strips path and extensition information from the
            # identifier to make an easily readable name.
            name = os.path.split(identifier)[-1]
            name = os.path.splitext(name)[0]

            for sub in self.subclasses:
                # Get the file for the identifier
                filepath = get_or_fetch(identifier, extensions='pdb.gz',
                                        urls=sub.urls)

                # If the file couldn't be found, raise an error and exit
                if filepath is None:
                    check_file(identifier, critical=True)

                # Load the file into a stream
                if filepath.endswith('.gz'):
                    with io.BufferedReader(gzip.open(filepath)) as stream:
                        content = stream.read().decode('latin-1')
                else:
                    with open(filepath) as f:
                        content = f.read()

                returned_molecules = sub.parse(content, name,
                                               source_molecules,
                                               *args, **kwargs)

                # If molecules were returned then it was successfully
                # parsed. No more parsing needed.
                if returned_molecules:
                    molecules += returned_molecules
                    break

        return molecules

    def parse(self, stream, name=None, source_molecules=None, *args, **kwargs):
        """Parse a data stream.
        
        Parameters
        ----------
        stream: io.Stream
            The data stream to parse.
        name: str, optional
            The name of the molecules to parse.
        source_molecule: molecule (:obj:`mollib.Molecule`), optional
            If specified, the molecule will have their data replaced with the 
            parsed data. Only the first model will be loaded.
        args: tuple, optional
            Positional arguments
        kwargs: dict, optional
            Keywork arguments
        
        Returns
        -------
        molecules: list of :obj:`mollib.Molecule` objects
        """
        return None
