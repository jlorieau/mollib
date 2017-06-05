"""
Interactions labels (str) are labels used to define the interaction between
atoms or groups of atoms. These utility functions are used to change between
interaction labels, label keys and atom groups.

Interaction labels (str) have the following features:

- Each group of atoms is identified by a subunit (optional), a residue
  number and an atom name.
  
  - ex: 'A.13CA' is the CA atom of the 13th residue in subunit A.

- If the subunit is not specified, the subunit 'A' is assumed.

- Atom interactions are related to each other through a delimiter
  (default: '-')
  
  - ex: '14N-H' corresponds to the interaction between the N and H atoms
    of residue 14 in subunit A.

- The basic '+1' and '-1' relative operators are supported to refer to next
  and previous residues.
  
  - ex: '13N-C-1' corresponds to the interaction between the N atom of
    residue 13 and the C atom of residue 12. Alternative, '13C-12N' could
    also be used.

- Wildcards are supported for atom names (default: '#')

  - ex: '13CA-HA#' corresponds to the CA atom of residue 13 and the HA2
    and HA3 atoms of residue 13.

Interaction types (str) have the following features:

- They correspond to a series of atom types without a stereospecific
  assignment. ex: The interaction type for '14N-H' and '14H-N' is 'N-H'.
  The interaction type for '35CA-HA', '14CA-HA2' and 'A.2HA-CA' is 'CA-HA'.

- They only use relative residue numbers. ex: The interaction type of
  '35N-34C' is 'N-C-1'
"""
# TODO: Add atom groups. ex: 35CA-(CG1CG2)
# TODO: Move to an nmr module.

import re
from itertools import chain, product, permutations
import logging

from .ordered_set import OrderedSet


# Once an interaction label is split, this regex parses it into its atom
# components, like subunit, residue number and atom name.
_re_label_str = (r'(?P<subunit>[A-Z]+\.)?'
                 r'(?P<number>\d*)'
                 r'(?P<name>[A-Z0-9]+)'
                 r'(?P<rel>[\-\+]\d+)?'
                 r'(?P<wildcard>.?)')
_re_label = re.compile(_re_label_str)


#: The following is a regex to match interaction labels
_re_interaction_str = (r'(?P<interaction>([A-Z]+\.)?'
                       r'(\d+[A-Z]+\d*\#?)'
                       r'(\-([A-Z]+\.)?[A-Z\d\-\+]+\#?)*)')


def _strip_str(string, strip_chars):
    """Strip the strip_chars from the given string."""
    return ''.join([j for j in string if j not in strip_chars])


def split_interaction_label(label, delimiter='-'):
    """Splits an interaction label into its component atom labels.


    .. note:: This function is designed to use relative atom references, like
              C-1 or N+1. These are used to reference atoms in previous or
              subsequent residues, respectively.

    Parameters
    ----------
    label: str
        The string identifier for the interaction. ex: A.14N-H
    delimiter: str
        The delimiter used for separating atom references.

    Returns
    -------
    pieces: list of str
        A list of atom label pieces

    Examples
    --------
    >>> split_interaction_label('14N-H')
    ['14N', 'H']
    >>> split_interaction_label('A.14H-13C')
    ['A.14H', '13C']
    >>> split_interaction_label('A.15N-16C-17H')
    ['A.15N', '16C', '17H']
    >>> split_interaction_label('A.14H-13C-N-1')
    ['A.14H', '13C', 'N-1']
    >>> split_interaction_label('A.14N-C-1')
    ['A.14N', 'C-1']
    >>> split_interaction_label('A.N+1-14C')
    ['A.N+1', '14C']
    """
    # Break the label into pieces
    pieces = label.split(delimiter)

    # Join the string piece for which a piece is a number. (ex: 'C-1' becomes
    # ['C', '1'] and should be joined back to 'C-1'
    l = []
    for count, item in enumerate(pieces):
        if item.isdigit() and count > 0:
            l[count - 1] = l[count - 1] + delimiter + item
        else:
            l.append(item)
    return l


def sort_func(label):
    """Generate a sort key for the given interaction_label.

    Parameters
    ----------
    label: str
        The string identifier for the interaction. ex: A.14N-H

    Returns
    -------
    tuple
        A tuple to be used to sort the interaction labels.

    Examples
    --------
    >>> sort_func('14N-H')
    ('N-H', 'A', 14, '14N-H')
    >>> sort_func('13C')
    ('C', 'A', 13, '13C')
    >>> sort_func('B.14N-13C')
    ('C-N', 'B', 13, 'B.14N-13C')
    >>> sort_func('B.35CA-HA2')
    ('CA-HA', 'B', 35, 'B.35CA-HA2')
    >>> sort_func('13CA-HA2')
    ('CA-HA', 'A', 13, '13CA-HA2')
    >>> sort_func('13CA-HA3')
    ('CA-HA', 'A', 13, '13CA-HA3')
    """
    # Generate a key from the label
    key = interaction_key(label)

    # The interaction type is from the atom names.  Strip numbers and '#' from
    # the atom names.
    atom_names = [_strip_str(i[2], '0123456789#') for i in key]
    interaction_type = '-'.join(atom_names)

    # Use the first atom's chain id residue number as the sort key.
    chain_id = key[0][0]
    res_number = key[0][1]

    return interaction_type, chain_id, res_number, label


def interaction_type(label, delimiter='-'):
    """Generate the string for the interaction type of a given
    interaction_label.

    .. note:: interaction types use relative residue numbers (ex: 'N-C-1')

    Parameters
    ----------
    label: str
        The string identifier for the interaction. ex: A.14N-H
    delimiter: str
        The delimiter used for separating atom references.

    Returns
    -------
    interaction_type
        A string for the interaction type.
    delimiter: str
        The delimiter used for separating atom references.

    Examples
    --------
    >>> interaction_type('14N-H')
    'N-H'
    >>> interaction_type('13C')
    'C'
    >>> interaction_type('B.35CA-HA2')
    'CA-HA'
    >>> interaction_type('35N-C-1')
    'N-C-1'
    >>> interaction_type('35N-34C')
    'N-C-1'
    >>> interaction_type('34C-35N')
    'C-N+1'
    """
    pieces = split_interaction_label(label, delimiter)

    # Pieces will have numbers and '#' characters in them. as long as these
    # don't end in '-1' or '+2', then the numbers and pounts can be stripped.
    processed = []
    last_residue = None
    for piece in pieces:
        match = _re_label.match(piece)
        if not match:
            continue
        d = match.groupdict()

        # Get the name and the relative identifier. Numbers should be stripped
        # only from the name
        name = _strip_str(d['name'], '0123456789#')
        number = int(d['number']) if d['number'] else None

        # See if there is a relative identifier
        # Pull out the relative identifier, like '+1' or '-2'
        rel = d['rel'] if isinstance(d['rel'], str) else ''

        # Determine whether a relative residue number should be replaced to the
        # atom name
        if (last_residue is not None and number is not None and
            last_residue != number):
            relative_number = number - last_residue
            rel = "{:+}".format(relative_number)
        last_residue = number

        processed.append(name + rel)

    return delimiter.join(processed)


def _key_sorting_function(i):
    """Function used to sort interaction keys."""
    # The 'H' character is renamed to 'h' in the sort function to guarantee
    # that protons are listed last in an interaction key.
    name = i[2] if i[2][0] != 'H' else 'h' + i[2][1:]
    return i[0], i[1], name


def interaction_key(label, default_subunit='A',
                    wildcard_char='#', wildtype_numbers=('1', '2', '3'),
                    sort=True):
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
    sort: bool, optional
        If True, the key results will be sorted. This ensures that the order
        of items are predictable.

    Returns
    -------
    tuple
        The tuple key for the tensor. ex: (('A', 14, 'N'), ('A', 14, 'H'))

    Examples
    --------
    >>> interaction_key('14N-H')
    (('A', 14, 'N'), ('A', 14, 'H'))
    >>> interaction_key('14N--H')  # Extra dashes are ignored
    (('A', 14, 'N'), ('A', 14, 'H'))
    >>> interaction_key('14H-N')  # The sort function orders the results
    (('A', 14, 'N'), ('A', 14, 'H'))
    >>> interaction_key('A.14H-13C')
    (('A', 13, 'C'), ('A', 14, 'H'))
    >>> interaction_key('B.18C')
    (('B', 18, 'C'),)
    >>> interaction_key('C')
    (('A', None, 'C'),)
    >>> interaction_key('13CA-HA#')
    (('A', 13, 'CA'), ('A', 13, 'HA1', 'HA2', 'HA3'))
    >>> interaction_key('35HB2')
    (('A', 35, 'HB2'),)
    >>> interaction_key('13N-C-1')
    (('A', 12, 'C'), ('A', 13, 'N'))
    >>> interaction_key('13C-N+1')
    (('A', 13, 'C'), ('A', 14, 'N'))
    """
    # Split the string into atom labels
    pieces = split_interaction_label(label)
    key = []

    # Parse each piece in terms of the residue number and atom name
    prev_residue_number = None
    prev_subunit = default_subunit

    for piece in pieces:
        # Match the string piece
        match = _re_label.match(piece)

        if not match:
            continue

        # Get the residue number and atom name
        groupdict = match.groupdict()
        residue_number = groupdict['number']
        atom_name = groupdict['name']
        subunit = groupdict['subunit']
        subunit = subunit.strip('.') if subunit else prev_subunit
        rel = groupdict['rel']
        wildcard = groupdict['wildcard']

        # If the residue number has not been specified (i.e. it is ''),
        # then determine the residue number from the previous piece
        if not residue_number:
            # use the previous residue number, unless a relative (rel) label
            # is specified like '-1' or '+1'
            if rel:
                msg = "Unable to parse relative residue identified in {}"
                try:
                    rel_int = int(rel)
                except ValueError:
                    logging.info(msg.format(label))
                    continue
                if prev_residue_number is None:
                    logging.info(msg.format(label))
                    continue
                residue_number = prev_residue_number + rel_int
            else:
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

    # If specified, sort the results. Either way, then turn the results into a
    # tuple to return
    if sort:

        key = sorted(key, key=_key_sorting_function)
        key = tuple(key)
    else:
        key = tuple(key)

    return key


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
        count = len([x for x in atom_names if x.startswith(group)])

        # Only process further if there are two or more atom_names that
        # match the group
        if count < 2:
            continue

        # Convert the atom_names into groupings, using an ordered_set
        l = [(i if not i.startswith(group)
              else _strip_str(i, '0123456789') + wildcard_char)
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

    .. note:: By default, this function will resort the order of items in the
              label to guarantee a predictable order.

    Parameters
    ----------
    label: str
        An interaction label. ex: 15N or 14N-H.
    key: bool, optional
        If True, the key results will be sorted. This ensures that the order
        of items are predictable.

    Returns
    -------
    label: str or None
        A reformatted, unambiguous interaction label. ex: A.15N or 14N-H.
        None is returned if the label cannot be parsed.

    Examples
    --------
    >>> validate_label('14N-H')
    'A.14N-H'
    >>> validate_label('A.14H-13C')
    'A.13C-14H'
    >>> validate_label('A.14H-B.13C')
    'A.14H-B.13C'
    >>> validate_label('B.18C')
    'B.18C'
    >>> validate_label('13CA-HA#')
    'A.13CA-HA#'
    >>> validate_label('35HB2')
    'A.35HB2'
    >>> validate_label('C')
    >>> validate_label('13N-C-1')
    'A.12C-13N'
    """
    try:
        key = interaction_key(label, *args, **kwargs)
        reformated_label = interaction_label(key, *args, **kwargs)
    except:
        msg = "The label '{}' does not correspond to an interaction label."
        logging.error(msg.format(label))
        return None
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


def interaction_atoms(key_or_label, molecule, sort=True):
    """A list of atom combinations for atoms identified by the interaction key
    or label.

    Parameters
    ----------
    molecule: :obj:`mollib.Molecule`
        A molecule from which combinations of atoms will be returned.
    key_or_label: str or tuple
        An interaction_label or interaction_key identifying atoms whose
        combinations are returned.
    sort: bool, optional
        If True, the key results will be sorted. This ensures that the order
        of items are predictable.

    Returns
    -------
    list of list of interaction atoms
        A list of combinations of interacting atoms.

    Examples
    --------
    >>> from mollib import Molecule
    >>> mol = Molecule('2MJB')
    >>> interaction_atoms('13C', mol)
    [[A.I13.C]]
    >>> interaction_atoms('35HA#', mol)
    [[A.G35.HA2], [A.G35.HA3]]
    >>> interaction_atoms('14N-H', mol)
    [[A.T14.N, A.T14.H]]
    >>> interaction_atoms('A.14H-13C', mol)
    [[A.I13.C, A.T14.H]]
    >>> interaction_atoms('35CA-HA#', mol)
    [[A.G35.CA, A.G35.HA2], [A.G35.CA, A.G35.HA3]]
    >>> interaction_atoms('33CA-HA-HB#', mol)
    [[A.K33.CA, A.K33.HA, A.K33.HB2], [A.K33.CA, A.K33.HA, A.K33.HB3]]
    >>> interaction_atoms('14N-C-1', mol)
    [[A.I13.C, A.T14.N]]
    """
    if isinstance(key_or_label, str) or isinstance(key_or_label, unicode):
        key = interaction_key(key_or_label, sort=sort)
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
    return [y for y in [[_key_to_atom(molecule, *x) for x in i]
                        for i in combs] if all(y)]


def get_dict_value(dictionary, label, delimiter='-', *args, **kwargs):
    """Return the value from the dictionary (grouped by interaction type) that
    matches the given interaction label string.

    This function is needed because the interaction label may be ordered
    differently than the interaction type in a dictionary. For example, the
    'A.14N-H' and 'A.14H-N' labels are both of the 'N-H' interaction type.

    Parameters
    ----------
    dictionary: dict
        The dictionary to get the value from.
        
        - **key**: The interaction type string. (ex: 'N-H' or 'CA-HA')
        - **value**: The value to return
    label: str
        The string identifier for the interaction. ex: A.14N-H
    delimiter: str
        The delimeter used for separating atom references.
    args: tuple, optional
        If specified a default argument, then this will be returned if the key
        isn't found. Otherwise a ValueError exception is raised.
    kwargs: dict, optional
        If specified a default argument, then this will be returned if the key
        isn't found. Otherwise a ValueError exception is raised.

    Raises
    ------
    ValueError
        If a default is not specified, then a ValueError is raised if the
        interaction type key is not found in the dictionary.

    Examples
    --------
    >>> get_dict_value({'N-H': 13.5}, '14N-H')
    13.5
    >>> get_dict_value({'N-H': 13.5}, '14H-N')
    13.5
    >>> get_dict_value({'CA-HA': 4.5}, '14N-H')
    Traceback (most recent call last):
    ...
    ValueError: Could not find the interaction 'N-H' in the given dictionary
    >>> get_dict_value({'CA-HA': 4.5}, '14N-H', default='default')
    'default'
    >>> get_dict_value({'N-C-1': 3}, '14N-C-1')
    3
    """
    # Get the interaction type.
    int_type = interaction_type(label, delimiter)  # ex: '14N-H' -> 'N-H'

    # We use the split_interaction label here so that relative atom numbers,
    # like 'N-C-1', aren't split. ex: 'N-C-1' -> ['N', 'C-1']
    pieces = split_interaction_label(int_type, delimiter)

    # Get all of the permutations of interaction types.
    # ex: 'N-H' -> ['N-H', 'H-N']
    key_permutations = [delimiter.join(i) for i in permutations(pieces)]

    # See if one of the permutations matches a key in the dict.
    for key in key_permutations:
        if key in dictionary:
            return dictionary[key]

    # At this point, the key wasn't found, so either return a default value or
    # raise an exception
    if len(args) > 0:
        return args[0]
    elif 'default' in kwargs:
        return kwargs['default']
    else:
        msg = "Could not find the interaction '{}' in the given dictionary"
        raise ValueError(msg.format(int_type))





