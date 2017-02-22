"""
These are utility functions for changing between string labels and label keys.
String labels are useful for saving text files, while label keys are useful
for identifying atoms, groups or atoms or other objects.
"""
import re
from string import digits
from itertools import chain
from .ordered_set import OrderedSet


re_label = re.compile(r'(?P<subunit>[A-Z]+\.)?'
                      r'(?P<number>\d*)'
                      r'(?P<name>[A-Z]+)'
                      r'(?P<wildcard>.?)')


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
    ('B', 18, 'C')
    >>> interaction_key('C')
    ('A', None, 'C')
    >>> interaction_key('13CA-HA#')
    (('A', 13, 'CA'), ('A', 13, 'HA1', 'HA2', 'HA3'))
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


    if len(key) == 1:
        return tuple(key[0])
    else:
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
    >>> interaction_label(('A', 13, 'C'))
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
    # Determine whether we have a multi key (i.e. a key that is a tuple of
    # tuples) or just a key (i.e. a simple tuple of items)
    multi_key = True if isinstance(key[0], tuple) else False

    # Check that None is not in the key
    if None in key or (multi_key and None in chain(*key)):
        raise KeyError('None cannot be in the key {}'.format(key))

    # If it's not a multi key, then simply format the label
    if not multi_key:
        subunit = key[0] + '.'
        res_number = key[1]
        atom_names = key[2:]
        return _convert_sublabel(wildcard_char, group_startswith, subunit,
                                 res_number, *atom_names)

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


def interaction_permutations(molecule, key_or_label):
    pass
