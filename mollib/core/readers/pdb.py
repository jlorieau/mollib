"""
Readers for Protein Databank Files.
"""
import re
from collections import deque
import itertools

import numpy as np

from .. import settings as settings
from .base import MoleculeReader
from mollib.utils.iteration import wrapit


def parse_atom_lines(atom_list, atom_cls):
    """Parse a list of 'ATOM' or 'HETATM' strings.

    Replaces the given list of strings with a list of atoms.

    Parameters
    ----------
    atom_list: list of str
        A list of 'ATOM' strings.
    atom_cls:
        The class to use in creating new atom objects.

    Returns
    -------
    A list of atom objects (:obj:`mollib.Atom`) if the 'ATOM' strings can be
    parsed.
    """
    processed = [atom_cls(number=int(l[6:11]),
                          name=str(l[12:16]).strip(),
                          residue=(str(l[17:20]).strip(),
                                   int(l[22:26]),),  # ('THR', 23)
                          chain=str(l[21:22]).strip() or 'A',
                          pos=np.array((float(l[31:38]),
                                        float(l[39:46]),
                                        float(l[47:54]),)),
                          element=str(l[76:78]).strip(),)
                 for l in atom_list]
    return processed


# The regex to match 'MODEL' lines in a pdb file.
_re_model = re.compile(r'^MODEL\s+(?P<model_id>[\d]+)\s+$')


def parse_model_line(model_line):
    """Parse a 'MODEL' line and return the model id number.
    
    Parameters
    ----------
    model_line: str
        The string that contains the new model information.
        ex: 'MODEL         1'
    
    Returns
    -------
    model_id: int,
        The parsed model_id number.
    """
    global _re_model

    match = _re_model.match(model_line)
    if match:
        model_id = int(match.groupdict()['model_id'])
        return model_id
    else:
        return None

# The regex to match 'CONECT' lines in a pdb file for bonding HETATM atoms.
_re_conect = re.compile(r'^CONECT(?P<numbers>[\s\d]+)$')


def parse_conect_lines(conect_lines):
    """Parser the 'CONECT' lines and return a connections list.
    
    Parameters
    ----------
    conect_lines: list of str
        The 'CONECT' lines of the PDB file.
        ex: CONECT 2059 2044
    
    Returns
    -------
    connections: list of tuple
        A list of tuples of atom numbers that should be connected.
    """
    global _re_conect

    connections = []
    for match in map(_re_conect.match, conect_lines):
        number_str = match.groupdict()['numbers']

        str_len = len(number_str)
        offset = -7  # This is for the stripped 'CONECT' string
        cols = [(7, 11), (12, 16), (17, 21), (22, 26), (27, 31), (32, 36),
                (37, 41), (42, 46), (47, 51), (52, 56), (57, 61)]
        cols = [(i + offset, j + offset + 1) for i, j in cols]

        atom_numbers = [number_str[i:j].strip()
                        if (i < str_len and j < str_len) else None
                        for i, j in cols]
        atom_numbers = [int(i) if i else None for i in atom_numbers]
        connections.append(atom_numbers)

    return connections


class PDBReader(MoleculeReader):
    """A molecule reader for protein databank (PDB) files.
    """

    order = 0

    def __init__(self, urls=None, *args, **kwargs):
        if urls is None:
            self.urls = settings.pdb_urls
        super(PDBReader, self).__init__(*args, **kwargs)

    def new_molecule(self, molecules, molecule_name, source_molecules=None,
                     model_id=None, model_ids=None):
        """Either retrieve of create a new molecule object.
        
        Parameters
        ----------
        molecules: list of molecule objects (:obj:`mollib.Molecule`)
            A list of the molecules that are being parsed and processed.
        molecule_name: str
            The name of the new molecule to create.
        source_molecules: `collections.deque` or :obj:`mollib.Molecule`
            If source_molecules are specified, these will be used instead of
            creating new molecules. Once the source_molecules run out, then None
            is returned.
        model_id: int, optional
            If specified, the new molecule will only be created if the model_id
            is in the model_ids specified. This will also assign the 'model_id'
            attribute for the molecule.
        model_ids: set of int, optional
            If specified a molecule will only be returned if the model_id is in
            the model_ids set.
        
        Returns
        -------
        molecule or None
            If a molecule is available or can be created, it will be returned.
            If a molecule is unavailable or cannot be created, None will be
            returned.
        """

        # First see if a model_id was specified and whether its in the allowed
        # model_ids. If not, return no molecules
        if (model_id is not None and
            model_ids is not None and
            model_id not in model_ids):
            return None

        # Now we have a valid model or molecule. Either we try to retrieve the
        # molecule from the source molecules or we create a new one.

        # First see if there are source_molecules to retrieve. If so, get one
        # of these. If they're depleted, then no more molecules can be returned
        if source_molecules is not None:
                if source_molecules:
                    # There are source molecules to return
                    molecule = source_molecules.popleft()
                    molecule.clear()
                    molecule.model_id = model_id
                    molecules.append(molecule)
                    return molecule
                else:
                    # Source molecules are depleted
                    return None

        # In this case, there are no source molecules and the molecule has an
        # allowed model_id (or it's not a model at all). Create a new molecule.
        molecule = self.molecule_class(identifier=molecule_name,
                                       use_reader=False)
        molecule.model_id = model_id
        molecules.append(molecule)
        return molecule

    def parse(self, stream, name=None, source_molecules=None,
              model_ids=None, *args, **kwargs):
        """Parse the file or string stream.

        Parameters
        ----------
        stream: stream
            The stream of text to parse.
        name: str, optional
            If specified, this is the name that will be given to the molecules
            and assigned to the molecule's 'name' attribute.
        source_molecules: `collections.deque` or :obj:`mollib.Molecule`
            If source_molecules are specified, these will be used instead of
            creating new molecules. Once the source_molecules run out, then None
            is returned.
        model_ids: list of int
            If specified, only the given model id numbers will be returned.
        args: tuple
            Function arguments
        kwargs: dict
            Function keywork arguments

        Returns
        -------
        molecules: list of molecule objects (:obj:`mollib.Molecule`)
            The molecules that we parsed and read in.
        """
        # Parse the arguments.
        molecule_name = name if isinstance(name, str) else 'Unnamed'

        source_molecules = (deque((source_molecules,))
                            if source_molecules is not None else None)

        # model_ids must be a list of integers
        model_ids = wrapit(model_ids)
        if model_ids is not None and not all(isinstance(i, int)
                                             for i in model_ids):
            msg = ("The specified model ids must be integer numbers. The "
                   "values '{}' were received.")
            raise TypeError(msg.format(model_ids))

        # Prepare the list of returned molecules. The keys are the model_id
        # numbers and the values are the molecules.
        molecules = []

        # Seperate the lines based on their type
        atom_lines = [line for line in stream.splitlines()
                      if line.startswith('ATOM') or
                      line.startswith('HETATM') or
                      line.startswith('MODEL') or
                      line.startswith('CONECT')]

        # Group contiguous lines
        groups = itertools.groupby(atom_lines, key=lambda x: x[0:6])

        # Placeholder values for the following loop
        current_model_id = None
        current_molecule = None
        connections = None

        for k, lines in groups:
            # Process 'MODEL' lines
            # This part checks to see if a new molecule should be assigned
            # based on a new model
            if k.startswith('MODEL'):
                # See if this is a new model and try to create it if needed.
                model_id = parse_model_line(next(lines))
                current_model_id = model_id
                current_molecule = self.new_molecule(molecules=molecules,
                    molecule_name=name, source_molecules=source_molecules,
                    model_id=model_id, model_ids=model_ids)
                continue

            # If a new molecule was not assigned based of a 'MODEL' line, then
            # see if a new one should be created.
            if current_molecule is None:
                current_molecule = self.new_molecule(molecules=molecules,
                    molecule_name=name, source_molecules=source_molecules,
                    model_id=current_model_id, model_ids=model_ids)

            # Parse the 'ATOM' and 'HETATM' lines (if a molecule has been
            # assigned
            if (current_molecule is not None and
                (k.startswith('ATOM') or k.startswith('HETATM'))):
                # Convert the 'ATOM' or 'HETATOM' lines into atoms
                # These are all part of the same molecule.
                atoms = parse_atom_lines(lines, atom_cls=self.atom_class)

                # Create the residue and chains as needed
                for atom in atoms:
                    # Get the chain_id string and residue number and name
                    # These will be replaced by the actual chain and residue
                    # objects. Also HETATM chains have a '*' next to them.
                    chain_id = (atom.chain if k.startswith('ATOM')
                                else atom.chain + '*')
                    res_name, res_num = atom.residue

                    # Retrieve or create the chain object
                    if chain_id not in current_molecule:
                        chain = self.chain_class(identifier=chain_id)
                        chain.molecule = current_molecule
                        current_molecule[chain_id] = chain
                    chain = current_molecule[chain_id]

                    # Retrieve or create the chain object
                    if res_num not in chain:
                        residue = self.residue_class(name=res_name,
                                                     number=res_num)
                        residue.chain = chain
                        residue.molecule = current_molecule
                        chain[res_num] = residue
                    residue = chain[res_num]

                    # Configure the atom
                    residue[atom.name] = atom

                    atom.chain = chain
                    atom.residue = residue
                    atom.molecule = current_molecule
                continue

            # Parse 'CONECT' lines
            if k.startswith('CONECT'):
                connections = parse_conect_lines(lines)

        # The parsing is finished.

        # Copy the connections across all molecules
        if connections is not None:
            for molecule in molecules:
                molecule.connections = connections

        # Now go through all of the molecules, and conduct needed processing
        # functions.
        for molecule in molecules:
            # Add prev_residue and next_residue pointers in the molecules.
            molecule.link_residues()

            # Translate atom names, if non-standard atom names are used.
            # translate_atom_name(molecule)

            # Set the atom topologies from the 'CONECT' records
            # TODO: This test should skip CONECT items if there are more than
            #       99,999 atoms
            molecule.set_atom_topologies()

        return molecules

#
# def translate_atom_name(molecule):
#     """Translate atom names that are not in the standard PDB format.
#
#     This is needed to process files from Xplor-NIH, for example, which uses 'HN'
#     and 'HB1'/'HB2' atom names instead of the standard 'H' and 'HB2'/'HB3'.
#     """
#     for residue in molecule.residues:
#         if 'HN' in residue:
#             atom = residue.pop('HN')
#             atom.name = 'H'
#             residue['H'] = atom
#         # if 'HT1' in residue and 'HT2' in residue and 'HT3' in residue:
#
#         if 'HA1' in residue and 'HA2' in residue:
#             atom1 = residue.pop('HA1')
#             atom2 = residue.pop('HA2')
#             atom1.name = 'HA2'
#             atom2.name = 'HA3'
#             residue['HA2'] = atom1
#             residue['HA3'] = atom2
