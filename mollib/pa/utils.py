import re
from itertools import chain


re_label = re.compile(r'(?P<subunit>[A-Z]+\.)?(?P<number>\d*)(?P<name>[A-Z]+)')

# TODO: rename to interaction_label.     >>> interaction_key('13CA-HA#')
#    ('A', 13, 'CA', 'A', 13, 'HA#')
def tensor_str_to_key(string, default_subunit='A'):
    """Convert a tensor string into a tuple key.

    Parameters
    ----------
    string: str
        The string identifier for the tensor. ex: A.14N-H

    Returns
    -------
    tuple
        The tuple key for the tensor. ex: ((14, 'N'), (14, 'H'))

    Examples
    --------
    >>> tensor_str_to_key('14N-H')
    (('A', 14, 'N'), ('A', 14, 'H'))
    >>> tensor_str_to_key('14N--H')
    (('A', 14, 'N'), ('A', 14, 'H'))
    >>> tensor_str_to_key('A.14H-13C')
    (('A', 14, 'H'), ('A', 13, 'C'))
    >>> tensor_str_to_key('B.18C')
    ('B', 18, 'C')
    >>> tensor_str_to_key('C')
    ('A', None, 'C')
    """
    # Split the string about the '-' character
    pieces = string.split('-')
    key = []

    # Parse each piece in terms of the residue number and atom name
    prev_residue_number = None
    prev_subunit = default_subunit

    for piece in pieces:
        # Match the string piece
        match = re_label.match(piece)

        if not match:
            continue

        # Get the residue number and atome name
        groupdict = match.groupdict()
        residue_number = groupdict['number']
        atom_name = groupdict['name']
        subunit = groupdict['subunit']
        subunit = subunit.strip('.') if subunit else prev_subunit

        # If the residue number has not been specified (i.e. it is ''),
        # then use the previous residue number even if it's equal to None
        if not residue_number:
            residue_number = prev_residue_number
        else:
            # Otherwise convert it to an integer
            residue_number = int(residue_number)

        # Append to the key and prepare for next loop iteration
        key.append((subunit, residue_number, atom_name))

        prev_residue_number = residue_number
        prev_subunit = subunit


    if len(key) == 1:
        return tuple(key[0])
    else:
        return tuple(key)



def tensor_key_to_str(key):
    """Convert a tensor tuple key to a tensor string.

    Parameters
    ----------
    key: tuple
        A tupled with the subunit (str), residue number (int) and atom name
        (str).

    Returns
    -------
    str:
        The tensor identifier string.

    Examples
    --------
    >>> tensor_key_to_str(('A', 13, 'C'))
    'A.13C'
    >>> tensor_key_to_str((('A', 14, 'N'), ('A', 14, 'H')))
    'A.14N-H'
    >>> tensor_key_to_str((('A', 14, 'N'), ('B', 15, 'H')))
    'A.14N-B.15H'
    >>> tensor_key_to_str((None, 13, 'C'))
    Traceback (most recent call last):
    ...
    KeyError: "None cannot be in the key (None, 13, 'C')"
    >>> tensor_key_to_str((('A', None, 'N'), ('A', 14, 'H')))
    Traceback (most recent call last):
    ...
    KeyError: "None cannot be in the key (('A', None, 'N'), ('A', 14, 'H'))"
    """
    # Check that None is not in the key
    if None in chain(key) or (isinstance(key[0], tuple) and
                              None in chain(*key)):
        raise KeyError('None cannot be in the key {}'.format(key))


    # Check to see if it's not a tuple of tuples (ex: ('A', 13, 'C'))
    if not isinstance(key[0], tuple):
        return "{}.{}{}".format(key[0], key[1], key[2])

    # In this case, it is a tuple of tuples.
    prev_subunit = None
    prev_res_number = None
    string_list = []

    for subunit, res_number, atom_name in key:
        # Format the string. The string tests are to avoid repitition in the
        # produced string.
        string = "{}{}{}".format((subunit + '.'
                                  if subunit != prev_subunit else ''),
                                 (res_number
                                  if res_number != prev_res_number else ''),
                                 atom_name)
        string_list.append(string)

        # Prepare for next loop iteration
        prev_subunit = subunit
        prev_res_number = res_number

    return '-'.join(string_list)



