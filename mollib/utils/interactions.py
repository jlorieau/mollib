"""
These are utility functions for changing between interaction labels and label
keys. String labels are useful for saving text files, while label keys are
useful for identifying interaction relationships between atoms.
"""
import re
from string import digits
from itertools import chain, product
from .ordered_set import OrderedSet


re_label = re.compile(r'(?P<subunit>[A-Z]+\.)?'
                      r'(?P<number>\d*)'
                      r'(?P<name>[A-Z0-9]+)'
                      r'(?P<wildcard>.?)')


class ValidationError(Exception): pass


def interaction_key(label, default_subunit='A',
                    wildcard_char='#', wildtype_numbers=('1', '2', '3')):
    """Convert an interaction label string into a tuple key.

    Keys either are simple keys with a tuple of (subunit, residue number,
    atom names), or they're multi-keys with tuple of tuples for these.

    Parameters
    ----------
    label: str
        The string identifier for the interaction. ex: A.14N-H
    default_subunit: str, optional
        If the subunit isn't specified in the label, use this subunit string
        as the default.
    wildcard_char: str, optional
        Whenever the wildcard_char is encountered in an atom name, it is
        interpreted as referring to multiple atoms. In this case, create an
        atom name in the key for each of the wildtype_numbers (below)
    wildtype_numbers: tuple, optional
        The atom names to create in the key whenever a wildtype_char is
        encountered.

    Returns
    -------
    tuple
        The tuple key for the tensor. ex: (('A', 14, 'N'), ('A', 14, 'H'))

    Examples
    --------
    >>> interaction_key('14N-H')
    (('A', 14, 'N'), ('A', 14, 'H'))
    >>> interaction_key('14N--H')
    (('A', 14, 'N'), ('A', 14, 'H'))
    >>> interaction_key('A.14H-13C')
    (('A', 14, 'H'), ('A', 13, 'C'))
    >>> interaction_key('B.18C')
    (('B', 18, 'C'),)
    >>> interaction_key('C')
    (('A', None, 'C'),)
    >>> interaction_key('13CA-HA#')
    (('A', 13, 'CA'), ('A', 13, 'HA1', 'HA2', 'HA3'))
    >>> interaction_key('35HB2')
    (('A', 35, 'HB2'),)
    """
    # Split the string about the '-' character
    pieces = label.split('-')
    key = []

    # Parse each piece in terms of the residue number and atom name
    prev_residue_number = None
    prev_subunit = default_subunit

    for piece in pieces:
        # Match the string piece
        match = re_label.match(piece)

        if not match:
            continue

        # Get the residue number and atom name
        groupdict = match.groupdict()
        residue_number = groupdict['number']
        atom_name = groupdict['name']
        subunit = groupdict['subunit']
        subunit = subunit.strip('.') if subunit else prev_subunit
        wildcard = groupdict['wildcard']

        # If the residue number has not been specified (i.e. it is ''),
        # then use the previous residue number even if it's equal to None
        if not residue_number:
            residue_number = prev_residue_number
        else:
            # Otherwise convert it to an integer
            residue_number = int(residue_number)

        # Append to the key and prepare for next loop iteration
        if wildcard == wildcard_char:
            items = [subunit, residue_number]
            items += [atom_name + str(i) for i in wildtype_numbers]
            key.append(tuple(items))
        else:
            key.append((subunit, residue_number, atom_name))

        prev_residue_number = residue_number
        prev_subunit = subunit

    return tuple(key)


def _convert_sublabel(wildcard_char, group_startswith,
                      subunit, residue_number, *atom_names):
    """(Private) helper function to convert a interaction key triplet
    into a string.

    Examples
    --------
    >>> _convert_sublabel('#', ('H',), 'A.', 25, 'HA1')
    'A.25HA1'
    >>> _convert_sublabel('#', ('H',), 'A.', 25, 'HA1', 'HA2')
    'A.25HA#'
    """
    # convert atom_names into a list
    atom_names = list(atom_names)

    for group in group_startswith:
        # Count the number of atom_names that start with this group of
        # startswith character
        count = len(filter(lambda x: x.startswith(group), atom_names))

        # Only process further if there are two or more atom_names that
        # match the group
        if count < 2:
            continue

        # Convert the atom_names into groupings, using an ordered_set
        l = [(i if not i.startswith(group)
              else i.translate(None, digits) + wildcard_char)
             for i in atom_names]
        atom_names_set = OrderedSet(l)
        atom_names = list(atom_names_set)

    # Now the atom_names have been compressed into wildcard_chars. Simply
    # format the sublabel
    return "{}{}".format(subunit, residue_number,) + ''.join(atom_names)


def interaction_label(key, wildcard_char='#', group_startswith=('H',)):
    """Convert an interaction key to a string.

    Keys either are simple keys with a tuple of (subunit, residue number,
    atom names), or they're multi-keys with tuple of tuples for these.

    Parameters
    ----------
    key: tuple
        A tupled with the subunit (str), residue number (int) and atom name
        (str).
    wildcard_char: str, optional
        Whenever multiple atom names match any of the group_startswith
        characters or string, group 2 or more of them with this wildcard_char.
    group_startswith: tuple, optional
        Atom names that begin with these characters will be grouped using
        the wildcard_char if there are 2 or more of them.

    Returns
    -------
    str:
        The tensor identifier label string.

    Examples
    --------
    >>> interaction_label((('A', 13, 'C'),))
    'A.13C'
    >>> interaction_label((('A', 14, 'N'), ('A', 14, 'H')))
    'A.14N-H'
    >>> interaction_label((('A', 14, 'N'), ('B', 15, 'H')))
    'A.14N-B.15H'
    >>> interaction_label((('A', 13, 'CA'), ('A', 13, 'HA1')))
    'A.13CA-HA1'
    >>> interaction_label((('A', 13, 'CA'), ('A', 13, 'HA1', 'HA2', 'HA3')))
    'A.13CA-HA#'
    >>> interaction_label((('A', None, 'N'), ('A', 14, 'H')))
    Traceback (most recent call last):
    ...
    KeyError: "None cannot be in the key (('A', None, 'N'), ('A', 14, 'H'))"
    """
    # Check that None is not in the key
    # if None in key or (multi_key and None in chain(*key)):
    if None in chain(*key):
        raise KeyError('None cannot be in the key {}'.format(key))

    # In this case, it's a multi-key that refers to multiple atoms.
    # Keep track of the previous subunit and residue number so that these
    # aren't inserted multiple times redundantly
    prev_subunit = None
    prev_res_number = None
    string_list = []

    for item in key:
        # Pull out the subunit and format it into a string. The test is to
        # remove redundant subunit characters in the formatted label
        subunit = item[0] + '.'
        subunit = subunit if subunit != prev_subunit else ''

        # Pull out the residue number and format it into a string. The test is
        # to remove redundant residue numbers in the formatted label
        res_number = str(item[1])
        res_number = res_number if res_number != prev_res_number else ''

        atom_names = item[2:]
        string = _convert_sublabel(wildcard_char, group_startswith,
                                   subunit, res_number, *atom_names)

        string_list.append(string)

        # Prepare for next loop iteration
        prev_subunit = subunit
        prev_res_number = res_number

    return '-'.join(string_list)


def validate_label(label, *args, **kwargs):
    """Reformat an interaction label so that it is unambiguous.

    Parameters
    ----------
    label: str
        An interaction label. ex: 15N or 14N-H.

    Returns
    -------
    label: str
        A reformatted, unambiguous interaction label. ex: A.15N or 14N-H.

    Raises
    ------
    ValidationError
        If the label is poorly formatted or there is too little information
        to be able to validate it.

    Examples
    --------
    >>> validate_label('14N-H')
    'A.14N-H'
    >>> validate_label('A.14H-13C')
    'A.14H-13C'
    >>> validate_label('A.14H-B.13C')
    'A.14H-B.13C'
    >>> validate_label('B.18C')
    'B.18C'
    >>> validate_label('13CA-HA#')
    'A.13CA-HA#'
    >>> validate_label('35HB2')
    'A.35HB2'
    >>> validate_label('C')
    Traceback (most recent call last):
    ...
    ValidationError: The label 'C' does not correspond to an interaction label.
    """
    try:
        key = interaction_key(label, *args, **kwargs)
        reformated_label = interaction_label(key, *args, **kwargs)
    except:
        msg = "The label '{}' does not correspond to an interaction label."
        raise ValidationError(msg.format(label))
    return reformated_label


def _key_to_atom(molecule, subunit, residue_number, atom_name):
    """(Private) helper function to convert a key to an atom for the given
    molecule."""
    if (subunit in molecule and
        residue_number in molecule[subunit] and
        atom_name in molecule[subunit][residue_number]):
        return molecule[subunit][residue_number][atom_name]
    else:
        return None


def interaction_atoms(key_or_label, molecule):
    """A list of atom combinations for atoms identified by the interaction key
    or label.

    Parameters
    ----------
    molecule: :obj:`mollib.Molecule`
        A molecule from which combinations of atoms will be returned.
    key_or_label: str or tuple
        An interaction_label or interaction_key identifying atoms whose
        combinations are returned.

    Returns
    -------
    list of list of interaction atoms
        A list of combinations of interacting atoms.

    >>> from mollib import Molecule
    >>> mol = Molecule('2MJB')
    >>> interaction_atoms('13C', mol)
    [[A.I13-C]]
    >>> interaction_atoms('35HA#', mol)
    [[A.G35-HA2], [A.G35-HA3]]
    >>> interaction_atoms('14N-H', mol)
    [[A.T14-N, A.T14-H]]
    >>> interaction_atoms('A.14H-13C', mol)
    [[A.T14-H, A.I13-C]]
    >>> interaction_atoms('35CA-HA#', mol)
    [[A.G35-CA, A.G35-HA2], [A.G35-CA, A.G35-HA3]]
    >>> interaction_atoms('33CA-HA-HB#', mol)
    [[A.K33-CA, A.K33-HA, A.K33-HB2], [A.K33-CA, A.K33-HA, A.K33-HB3]]

    """
    if isinstance(key_or_label, str):
        key = interaction_key(key_or_label)
    elif isinstance(key_or_label, tuple):
        key = key_or_label
    else:
        msg = ("The key or label '{}' must either be an identifier key or label"
               " string")
        raise KeyError(msg.format(key_or_label))

    # Prepare the combinations
    item_lists = []
    for item in key:
        subunit = item[0]
        res_number = item[1]
        atom_names = item[2:]

        item_lists.append([(subunit, res_number, i) for i in atom_names])

    combs = product(*item_lists)

    # Convert keys to atom instances.
    return filter(lambda y: all(y),
                  [map(lambda x: _key_to_atom(molecule, *x), i)
                   for i in combs])