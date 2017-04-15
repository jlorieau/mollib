# -*- coding: utf-8 -*-
"""
Read in RDC and RACS data from files.
"""
import re
import logging
import gzip

from mollib.utils.interactions import validate_label, interaction_label
from .data_types import RDC, RACS


def read_pa_file(filename):
    """Read data from a partial alignment data file.

    Parameters
    ----------
    filename: str
        The filename to read the data from. The file can be a gzipped file.

    Returns
    -------
    data: dict
        A dict with the data. The keys are tensor keys
        and the values are :obj:`RDC` or :obj:`RACS` datum objects.
    """
    data = {}

    if filename.endswith('.gz'):
        with gzip.open(filename) as f:
            string = f.read()
    else:
        with open(filename) as f:
            string = f.read()

    retrieved_data = read_pa_string(string)
    data.update(retrieved_data)

    # If reading the file as a pa string didn't work, try the DC file
    # format instead
    if len(retrieved_data) == 0:
        retrieved_data = read_dc_string(string)
        data.update(retrieved_data)

    # If reading the file as a DC file format didn't work, try the mr file
    # format instead
    if len(retrieved_data) == 0:
        retrieved_data = read_mr_string(string)
        data.update(retrieved_data)

    return data


re_pa = re.compile(r'^\s*'
                   r'(?P<interaction>[\w\-\.]+#?)'
                   r'\s+'
                   r'(?P<value>[eE\d\-\+\.]+)'
                   r'\s*'
                   r'(?P<error>[eE\d\-\+\.]*)')


def read_pa_string(string):
    """Read data from a partial alignment data string.

    Parameters
    ----------
    string: str or list of str
        Either a (multiline) string or a list of strings.


    .. note::

        The format of the file is as follows:

        ::

            # Interaction   Value (Hz)   Error (optional)
            14N-H           -14.5        0.1
            15N-H             3.5
            A.16N-H          -8.5        0.2  # larger error

            A.16H-A.15C       0.5        0.1
            B.16H-B.15C       0.5        0.1

            # Residual anisotropic chemical shift data
            # Interaction   Value (ppb)   Error (optional)
            5C                112         1
            6C               -250


    Returns
    -------
    data: dict
        A dict with the data. 
        - **key**: interaction labels (str). ex: '14N-H'
        - **value**: RDC or RACS datum objects (:obj:`RDC` :obj:`RACS`)

    Examples
    --------
    >>> data = read_pa_string('''
    ... # Interaction   Value (Hz)   Error (optional)
    ... 14N-H           -14.5        0.1
    ... 15N-H             3.5
    ...
    ... 5C                112         1''')
    >>> for k, v in sorted(data.items()):
    ...     print("{:<10} {}".format(k, v))
    A.14N-H    rdc(-14.5±0.1)
    A.15N-H    rdc(3.5±0.0)
    A.5C       racs(112.0±1.0)
    """
    # Convert the string into a list of lines, if it isn't already.
    if not isinstance(string, list):
        string = string.splitlines()

    # Prepare the returned data list
    data = {}

    # Find all of the matches and produce a generator
    matches = (m for l in string for m in [re_pa.search(l)] if m)

    # iterate over the matches and pull out the data.
    for match in matches:
        d = match.groupdict()
        logging.debug("read_pa_string match: " + str(d))

        interaction_key = validate_label(d['interaction'])
        value = float(d['value'])
        error = float(d['error'] if d['error'] else 0.0)

        # Add the Datum to the data list. If the interaction_label has a '-'
        # character, then it is referring the multiple atoms and must be a
        # residual dipolar coupling (RDC). Otherwise, it's a residual
        # anisotropic chemical shift (RACS).
        hyphen_count = interaction_key.count('-')
        if hyphen_count == 1:
            data[interaction_key] = RDC(value=value, error=error)
        elif hyphen_count == 0:
            data[interaction_key] = RACS(value=value, error=error)
        else:
            continue

    return data


re_dc = re.compile(r'^\s*'
                   r'(?P<res_num1>\d+)'
                   r'\s+[A-Z]{3}\s+'
                   r'(?P<atom_name1>[A-Z0-9]+#?)'
                   r'\s+'
                   r'(?P<res_num2>\d+)'
                   r'\s+[A-Z]{3}\s+'
                   r'(?P<atom_name2>[A-Z0-9]+#?)'
                   r'\s+'
                   r'(?P<value>[eE\d\-\+\.]+)')


def read_dc_string(string):
    """Read data from a DC RDC data string.

    Parameters
    ----------
    string: str or list of str
        Either a (multiline) string or a list of strings.

    Returns
    -------
    data: dict
        A dict with the data. 
        - **key**: interaction labels (str). ex: '14N-H'
        - **value**: RDC  datum objects (:obj:`RDC`)
    """
    # Convert the string into a list of lines, if it isn't already.
    if not isinstance(string, list):
        string = string.splitlines()

    # Prepare the returned data list
    data = {}

    # Find all of the matches and produce a generator
    matches = (m for l in string for m in [re_dc.search(l)] if m)

    # iterate over the matches and pull out the data.
    for match in matches:
        d = match.groupdict()
        logging.debug("read_dc_string match: " + str(d))

        # Get the residue numbers and atom names
        res_num1 = int(d['res_num1'])
        res_num2 = int(d['res_num2'])
        atom_name1 = d['atom_name1']
        atom_name2 = d['atom_name2']
        value = float(d['value'])

        # Generation the interaction key
        key = (('A', res_num1, atom_name1),
               ('A', res_num2, atom_name2))

        # Get the interaction label
        label = interaction_label(key)
        label = validate_label(label)

        # Add it to the dict
        data[label] = RDC(value=value)

    return data


re_mr = re.compile(r'assign'
                   r'\s*'
                   r'\(\s*resid\s*(?P<coord_num>\d+)[\s\w]+OO\)'
                   r'(?:\n[\s\w]*\([\w\s]+[XYZ]\s*\))+'
                   r'(?:\n[\s\w]*\([a-zA-Z\s]+'
                       r'(?P<res_i>\d+)[a-z\s]+'
                       r'(?P<name_i>[A-Z0-9\#]+)\s*\))?'
                   r'\s*'
                   r'(?P<value_i>[\d\.\-eE]+)?'
                   r'(?:\n[\s\w]*\([a-zA-Z\s]+'
                       r'(?P<res_j>\d+)[a-z\s]+'
                       r'(?P<name_j>[A-Z0-9\#]+)\s*\))'
                   r'\s*'
                   r'(?P<value_j>[\d\.\-eE]+)?',
                   re.MULTILINE)


def read_mr_string(string, coordinate_num=None):
    """Read data from a mr data string.

    .. note:: The current version only supports the 'A' subunit.

    Parameters
    ----------
    string: str 
        A multiline string.

    Returns
    -------
    data: dict
        A dict with the data. 
        - **key**: interaction labels (str). ex: '14N-H'
        - **value**: RDC  datum objects (:obj:`RDC`)
    """
    # Prepare the returned data list
    data = {}
    data_keys = {}

    # Find all of the matches and produce a generator
    matches = re_mr.finditer(string)

    # iterate over the matches and pull out the data.
    for match in matches:
        d = match.groupdict()
        logging.debug("read_mr_string match: " + str(d))

        # Get the residue numbers and atom names
        coord_num = d['coord_num']
        res_i = int(d['res_i'])
        name_i = d['name_i']
        res_j = int(d['res_j'])
        name_j = d['name_j']
        value = float(d['value_j'])

        # Only collect data with the matching coordinate_num (if one has been
        # specified either in the function parameters or in a previous
        # iteration of the loop)
        if coord_num != coordinate_num:
            if coordinate_num is None:
                coordinate_num = coord_num
            else:
                continue

        # if res_i/name_i are None then this is a CSA value. Generate the key
        # and data type
        if res_i is None:
            key = (('A', res_j, name_j),)
            data_type = RACS
        else:
            key = (('A', res_i, name_i), ('A', res_j, name_j))
            data_type = RDC

        # Add it to the dict
        data_keys[key] = data_type(value=value)

    # Convert the data_keys, which is identified by key tuples, to labels. This
    # is done separately so that incorrect atom names, like 'HN', can be
    # converted to standard PDB format, i.e. 'H'.

    for key, value in data_keys.items():
        new_key = []
        # Check to see if any atom names are
        for subunit, res_num, atom_name in key:
            if atom_name == 'HN':
                atom_name = 'H'
            if atom_name == 'HA1':
                atom_name = 'HA3'
            new_key.append((subunit, res_num, atom_name))

        new_key = tuple(new_key)

        # Get the interaction label
        label = interaction_label(new_key)
        label = validate_label(label)

        data[label] = value

    return data


