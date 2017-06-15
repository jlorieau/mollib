"""
Readers for Protein Databank Files.
"""
import re
from collections import OrderedDict, deque
import itertools

import numpy as np

from .. import settings as settings
from .base import MoleculeReader
from mollib.utils.iteration import wrapit


### New Implementation

def parse_atom_lines(atom_list, atom_cls):
    """Parse a list of 'ATOM' or 'HETATM' strings.

    A list of atoms replaces the list of strings passed to the function.

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
                          chain=str(l[21:22]).strip(),
                          pos=np.array((float(l[31:38]),
                                        float(l[39:46]),
                                        float(l[47:54]),)),
                          element=str(l[76:78]).strip(),)
                 for l in atom_list]
    return processed


# The regex to match 'MODEL' lines in a pdb file.
_re_model = re.compile(r'^MODEL\s+(?P<model_id>[\d]+)\s+$')


def parse_model_line(model_line):
    """Parse a 'MODEL' line and return the model id number."""
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
    """Parser the 'CONECT' lines and return a connections list."""
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

    order = -10

    def __init__(self, urls=None, *args, **kwargs):
        if urls is None:
            self.urls = settings.pdb_urls
        super(PDBReader, self).__init__(*args, **kwargs)

    def new_molecule(self, molecules, molecule_name, source_molecules=None,
                     model_id=None, model_ids=None):
        """Either retrieve of create a new molecule object.
        
        Parameters
        ----------
        
        Returns
        -------
        molecule or None
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
        molecule = self.molecule_class(name=molecule_name)
        molecule.model_id = model_id
        molecules.append(molecule)
        return molecule

    def parse(self, stream, name=None, source_molecules=None,
              model_ids=None, *args, **kwargs):
        """Parse the file stream.

        Parameters
        ----------
        stream
        name
        model_ids: list of int
            If specified, only the given model id numbers will be returned.
        args
        kwargs

        Returns
        -------

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

        # If model_ids aren't specified, see if there are any settings that
        # might influence which models to select
        if model_ids is None and settings.pdb_first_model:
            model_ids = [1, ]

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

            # Set the atom topologies from the 'CONECT' records
            # TODO: This test should skip CONECT items if there are more than
            #       99,999 atoms
            molecule.set_atom_topologies()

        return molecules


### Old Implementation

# TODO: Make a faster reader with structs and binary, with a fallback
#       to the regex matcher.

class PDBRigidRegexReader(MoleculeReader):
    """A Molecule Reader for Protein Databank (PDB) files.
     
    PDBRigidRegexReader is a single-threaded Pythong implementation using a
    rigid structure regex that strictly follows the column width
    specifications of PDB files.
    """

    order = 0

    def __init__(self, urls=None, *args, **kwargs):
        if urls is None:
            self.urls = settings.pdb_urls
        super(PDBRigidRegexReader, self).__init__(*args, **kwargs)

    # The regex to match 'ATOM' and 'HETATM' lines. This regex is strict in
    # terms of the positioning of items.
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

    # The following are the conversion functions for the different matched items
    # in the _re_atom regex.
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

    # The regex to match 'CONECT' lines in a pdb file for bonding HETATM atoms.
    _re_conect = re.compile(r'^CONECT(?P<numbers>[\s\d]+)$')

    # The regex to match 'MODEL' lines in a pdb file.
    _re_model = re.compile(r'^MODEL\s+(?P<model_id>[\d]+)\s+$')

    def _match_atom(self, match, molecule):
        """Matches an ATOM or HETATM line in a PDB file and adds the matched
        atom to the molecule.

        Parameters
        ----------
        match: :obj:`re.MatchObject`.
            A regex match object generated by self._re_atom
        molecule: :obj:`mollib.Molecule`
            The molecule to add the matched atom to.

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
        if identifier not in molecule:
            chain = self.chain_class(identifier=identifier)
            chain.molecule = molecule
            molecule[identifier] = chain

        chain = molecule[identifier]

        # create Residue, if it doesn't already exist
        res_number, res_name, = t[6], t[4]

        if res_number not in chain:
            try:
                residue = self.residue_class(number=res_number, name=res_name)

                residue.chain = chain
                residue.molecule = molecule
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
                               molecule=molecule)

        residue[name] = atom

    def _match_conect(self, match, molecule):
        """Matches a 'CONECT' line in a PDB file and populates the connections
        attribute.

        Parameters
        ----------
        match: :obj:`re.MatchObject`.
            A regex match object generated by self._re_conect
        molecule: :obj:`mollib.Molecule`
            The molecule to add the matched atom to.
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
                        if (i < str_len and j < str_len) else None
                        for i,j in cols]
        atom_numbers = [int(i) if i else None for i in atom_numbers]

        # The connections attribute is created by the molecule's __init__
        # function
        molecule.connections.append(atom_numbers)

    def _match_model(self, match, molecule):
        """Matches a MODEL line in a PDB file and sets the molecule's 'model_id'
        attribute.

        Parameters
        ----------
        match: :obj:`re.MatchObject`.
            A regex match object generated by self._re_model
        molecule:
            The molecule to set the match_id to.
        
        Returns
        -------
        model_id: int
            The model_id of the matched model
        """
        model_id = int(match.groupdict()['model_id'])
        molecule.model_id = model_id
        return model_id

    # The matchers have the starting string of matched items as keys, as a
    # tuple of the compiled regex and match processing functions.
    _matchers = OrderedDict((('ATOM', (_re_atom, _match_atom)),
                             ('MODEL', (_re_model, _match_model)),
                             ('CONECT', (_re_conect, _match_conect)),
                             ))

    def parse(self, stream, name=None, source_molecules=None,
              model_ids=None, *args, **kwargs):
        """
        
        Parameters
        ----------
        stream
        name
        model_ids: list of int
            If specified, only the given model id numbers will be returned.
        args
        kwargs

        Returns
        -------

        """
        # Parse the arguments.
        molecule_name = name if isinstance(name, str) else 'Unnamed'

        source_molecules = (deque((source_molecules, ))
                            if source_molecules is not None else None)

        # model_ids must be a list of integers
        model_ids = wrapit(model_ids)
        if model_ids is not None and not all(isinstance(i, int)
                                             for i in model_ids):
            msg = ("The specified model ids must be integer numbers. The "
                   "values '{}' were received.")
            raise TypeError(msg.format(model_ids))

        # If model_ids aren't specified, see if there are any settings that
        # might influence which models to select
        if model_ids is None and settings.pdb_first_model:
            model_ids = [1, ]

        # A list of regex matchers to harvest data from each line.
        # The first item is the regex to match. If there is a match,
        # the second item (function) will be called with the match
        matchers = self._matchers.copy()

        # Prepare the list of returned molecules
        molecules = []
        current_molecule = None

        # Find the ATOM/HETATM lines and pull out the necessary data
        for line in stream:
            # Placeholders for the regex match and the match processing function
            m, func = None, None

            # Gzipped files return bytes lines that have to be decoded
            if type(line) == bytes:
                line = line.decode('latin-1')

            # Advance the regex matching iterator
            for name, (regex, func) in matchers.items():
                m = regex.match(line)

                # If there's no match, nothing else can be done
                if not m:
                    continue

                # At this point, there is a match.

                # Set the current_molecule, if a molecule hasn't already been
                # specified (current_molecule is None) or a new model has been
                # specified.
                if m and (name == 'MODEL' or current_molecule is None):
                    # Use source molecules, if available
                    # if source_molecules:
                    #     current_molecule = source_molecules.popleft()
                    # else:
                    #     current_molecule = self.molecule_class(molecule_name,
                    #                                            use_reader=False)
                    if source_molecules is None:
                        current_molecule = self.molecule_class(molecule_name,
                                                               use_reader=False)
                        molecules.append(current_molecule)
                    elif source_molecules:
                        current_molecule = source_molecules.popleft()
                        molecules.append(current_molecule)
                    else:
                        del matchers['ATOM']
                        del matchers['MODEL']
                break

            # If line cannot be matched, skip this line.
            if m is None or func is None:
                continue

            # Parse the matched line.
            func(self, m, current_molecule)

        # The parsing is finished.

        # Copy the connections across all molecules
        if molecules:
            connections = molecules[-1].connections
            for molecule in molecules:
                molecule.connections = connections

        # Select the molecules by model_id, if needed
        if model_ids is not None and source_molecules is None:
            molecules = [m for m in molecules
                         if hasattr(m, 'model_id') and m.model_id in model_ids]

        # Now go through all of the molecules, and conduct needed processing
        # functions.
        for molecule in molecules:
            # Add prev_residue and next_residue pointers in the molecules.
            molecule.link_residues()

            # Set the atom topologies from the 'CONECT' records
            # TODO: This test should skip CONECT items if there are more than
            #       99,999 atoms
            molecule.set_atom_topologies()

        return molecules

