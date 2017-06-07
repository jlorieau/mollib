# -*- coding: utf-8 -*-
"""
Read in RDC and RACS data from files.
"""
import re
import logging
import gzip
from collections import OrderedDict

from mollib.utils.interactions import (validate_label, interaction_label,
                                       _re_interaction_str)
from .data_types import RDC, RACS
from . import logs


def get_data(data_sets, set_id=None):
    """Return the specified data set from multiple sets of data.
    
    Parameters
    ----------
    data_sets: :obj:`collections.OrderedDict` or dict
        The data sets (one or multiple) to retrieve a single data from.
        If a single data set is passed, then it will be returned, unchanged.
    set_id: str (optional)
        If specified, the data set matching the given set_id will be returned.
        The set_id can either be a string of the set_id or a number for the
        set number to return, starting from '0'.
        If not specified, the first data set will be returned
    
    Returns
    -------
    data: dict or None
        The selected dataset, or None if not found.
        
    Examples
    --------
    >>> from collections import OrderedDict
    >>> data_sets = OrderedDict((('A', {'14N-H': 1}), ('B', {'14N-H': 2})))
    >>> get_data(data_sets)
    {'14N-H': 1}
    >>> get_data(data_sets, 'B')
    {'14N-H': 2}
    >>> get_data(data_sets, 1)
    {'14N-H': 2}
    >>> get_data(data_sets['B'])
    {'14N-H': 2}
    """
    # Directly retrieve the dataset, if specified
    if set_id is not None:
        if set_id in data_sets:
            return data_sets[set_id]

        # Try to see if the set_id can be converted to a number that can be
        # used to access the set
        if isinstance(set_id, int):
            pass
        elif ((isinstance(set_id, str) or isinstance(set_id, unicode)) and
              set_id.isdigit()):
            set_id = int(set_id)
        else:
            msg = "Could not parse data set id '{}'".format(set_id)
            if msg not in logs.errors:
                logging.error(msg)
                logs.errors.add(msg)
            return None

        # The set_id is not an integer. Let's see if it corresponds to a key
        # in the key array
        keys = list(data_sets.keys())

        if set_id < len(keys):
            key = keys[set_id]
            return data_sets[key]

        # The set_id could not be parsed into a useful value to retrieve the
        # data set. Return None. (i.e. not found)
        msg = "Could not find data set '{}'".format(set_id)
        if msg not in logs.errors:
            logging.error(msg)
            logs.errors.add(msg)
        return None

    # At this point, no set_id key was specified. Just return the first
    # dataset, if available.
    keys = list(data_sets.keys())
    if len(keys) == 0:
        return None
    key = keys[0]

    # Return the item if its a dict, (i.e. a data set), otherwise, just
    # return the data_sets itself.
    if isinstance(data_sets[key], dict):
        return data_sets[key]
    else:
        return data_sets


# Structures known to not work: 2LZF (CYANA), 2LT9 (CYANA)
def read_pa_file(filename, set_id=None, ignore_ext=False):
    """Read data from a partial alignment data file.

    Parameters
    ----------
    filename: str
        The filename to read the data from. The file can be a gzipped file.
    set_id: str (optional)
        If specified, the data set matching the given set_id will be returned.
        The set_id can either be a string of the set_id or a number for the
        set number to return, starting from '0'.
        If not specified, the first data set will be returned
    ignore_ext: bool (optional)
        If True, all the parsers will be attempted for the file. Otherwise,
        the only specific parsers will be attempted. So far, this only works
        with '.mr' and '.mr.gz' files.

    Returns
    -------
    data: dict
        A dict with the data. The keys are tensor keys
        and the values are :obj:`RDC` or :obj:`RACS` datum objects.
    """
    data = OrderedDict()

    # Read in the filename to a string
    if filename.endswith('.gz'):
        with gzip.open(filename) as f:
            string = f.read()
    else:
        with open(filename) as f:
            string = f.read()

    # Convert the string, if it's in bytes, to a text string. This is needed
    # for Python 3 compatibility.
    if type(string) == bytes:
        string = string.decode('latin-1')

    # The objective here is to read from the most specific to the least
    # specific. First start with '.mr' data format.
    retrieved_data = read_mr_string(string, set_id)

    # Update the data dict if data was found
    if retrieved_data is not None and len(retrieved_data) > 0:
        data.update(retrieved_data)
        return data

    # If it's a .mr or .mr.gz file and ignore_ext is False, then no other
    # parsers should be attempted
    if not ignore_ext and (filename.endswith('mr') or
                           filename.endswith('.mr.gz')):
        return data

    # Now, try the next parser. Try the .pa data format
    retrieved_data = read_pa_string(string)

    # Update the data dict if data was found
    if retrieved_data is not None and len(retrieved_data) > 0:
        data.update(retrieved_data)
        return data

    # Now try NMRPipe DC file data format. This format is promiscuous and may
    # match inadvertantly.
    retrieved_data = read_dc_string(string)

    # Update the data dict if data was found
    if retrieved_data is not None and len(retrieved_data) > 0:
        data.update(retrieved_data)
        return data

    # If nothing else, return the empty dataset
    return data


_re_pa = re.compile(r'^\s*' +
                    _re_interaction_str +
                   r'\s+'
                   r'(?P<value>[\d\-\+\.]+[eE]?[\d\-\+]*)'
                   r'\s*'
                   r'(?P<error>[\d\-\+\.]*[eE]?[\d\-\+]*)')


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
    """
    # Convert the string into a list of lines, if it isn't already.
    if not isinstance(string, list):
        string = string.splitlines()

    # Prepare the returned data list
    data = OrderedDict()

    # Find all of the matches and produce a generator
    matches = (m for l in string for m in [_re_pa.search(l)] if m)

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


_re_dc = re.compile(r'^\s*'
                   r'(?P<res_num1>\d+)'
                   r'\s+[A-Z]{3}\s+'
                   r'(?P<atom_name1>[A-Z0-9]+#?)'
                   r'\s+'
                   r'(?P<res_num2>\d+)'
                   r'\s+[A-Z]{3}\s+'
                   r'(?P<atom_name2>[A-Z0-9]+#?)'
                   r'\s+'
                   r'(?P<value>[\d\-\+\.]+[eE]?[\d\-\+]*)')


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
    data = OrderedDict()

    # Find all of the matches and produce a generator
    matches = (m for l in string for m in [_re_dc.search(l)] if m)

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

        # Generate the interaction key
        key = (('A', res_num1, atom_name1),
               ('A', res_num2, atom_name2))

        # Get the interaction label
        label = interaction_label(key)
        label = validate_label(label)

        # Add it to the dict
        data[label] = RDC(value=value)

    return data


_re_mr = re.compile(r'assign'
                   r'\s*'
                   r'\([\s\w]*?resid?\s*(?P<coord_num>\d+)[\s\w]+OO\s*\)\s*'
                   r'(?:\n[\s\w]*\([\w\s]+[XYZ]\s*\)\s*){3}'
                   r'(\n[\s\w]*\(\s*'
                       r'(?:segid?\s+(?P<chain_i>[A-Z])\s+)?[\s\w]*?'
                       r'resid?\s+(?P<res_i>\d+)[a-z\s]+?'
                       r'name\s+(?P<name_i>[A-Z0-9\#]+)\s*\))?'
                   r'\s*'
                   r'(?P<value_i>[\d\-\+\.]+[eE]?[\d\-\+]*)?'
                   r'(\n[\s\w]*\(\s*'
                       r'(?:segid?\s+(?P<chain_j>[A-Z])\s+)?[\s\w]*?'
                       r'resid?\s*(?P<res_j>\d+)[a-z\s]+?'
                       r'name\s+(?P<name_j>[A-Z0-9\#]+)\s*\))'
                   r'\s*'
                   r'(?P<value_j>[\d\-\+\.]+[eE]?[\d\-\+]*)?',
                    re.MULTILINE)


# TODO: fix reading of 5T1N
def read_mr_string(string, set_id=None):
    """Read data from a mr data string.

    .. note:: The current version only supports the 'A' subunit.

    Parameters
    ----------
    string: str 
        A multiline string.
    set_id: str (optional)
        See :func:`get_data`.

    Returns
    -------
    data: dict
        A dict with the data. 
        - **key**: interaction labels (str). ex: '14N-H'
        - **value**: RDC  datum objects (:obj:`RDC`)
    """
    # Prepare the returned data list
    data_sets = OrderedDict()

    # Find all of the matches and produce a generator
    matches = _re_mr.finditer(string)

    # iterate over the matches and pull out the data.
    for match in matches:
        d = match.groupdict()
        logging.debug("read_mr_string match: " + str(d))

        # Get the residue numbers and atom names
        coord_num = d['coord_num']
        chain_i = d['chain_i'] if d['chain_i'] is not None else 'A'
        res_i = int(d['res_i']) if d['res_i'] is not None else None
        name_i = d['name_i']

        chain_j = d['chain_j'] if d['chain_j'] is not None else 'A'
        res_j = int(d['res_j']) if d['res_j'] is not None else None
        name_j = d['name_j']
        value = float(d['value_j']) if d['value_j'] is not None else None

        # Either one or the other residue must be specified as well as a value
        if (res_i is None and res_j is None) or value is None:
            continue

        # if res_i/name_i are None then this is a CSA value. Generate the key
        # and data type
        if res_i is None:
            key = ((chain_j, res_j, name_j),)
            data_type = RACS
        else:
            key = ((chain_i, res_i, name_i), (chain_j, res_j, name_j))
            data_type = RDC

        # Convert the key, which is a tuple into and interaction label. This
        # is done separately so that incorrect atom names, like 'HN', can be
        # converted to standard PDB format, i.e. 'H'.
        new_key = []
        for subunit, res_num, atom_name in key:
            if atom_name == 'HN':
                atom_name = 'H'
            if atom_name == 'HA1':
                atom_name = 'HA3'
            new_key.append((subunit, res_num, atom_name))
        new_key = tuple(new_key)
        label = interaction_label(new_key)
        label = validate_label(label)

        # Add it to the dict
        data = data_sets.setdefault(coord_num, OrderedDict())
        data[label] = data_type(value=value)

    # Return only the data set specified
    return get_data(data_sets, set_id)
