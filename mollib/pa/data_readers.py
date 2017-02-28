# -*- coding: utf-8 -*-
"""
Read in RDC and RACS data from files.
"""
import re

from mollib.utils.interactions import validate_label
from .data_types import RDC, RACS


re_pa = re.compile(r'\s*'
                   r'(?P<interaction>[\w\-\.]+)'
                   r'\s+'
                   r'(?P<value>[E\d\-\+\.]+)'
                   r'\s*'
                   r'(?P<error>[E\d\-\+\.]*)')


def read_pa_data(filename):
    """Read data from a partial alignment data file.

    Parameters
    ----------
    filename: str
        The filename to read the data from.

    Returns
    -------
    data: dict
        A dict with the data. The keys are tensor keys
        and the values are :obj:`RDC` or :obj:`RACS` datum objects.
    """
    data = {}
    with open(filename, 'r') as f:
        lines = list(f.readlines())
        data.update(read_pa_string(lines))
    return data


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
        A dict with the data. The keys are interaction keys
        and the values are :obj:`RDC` or :obj:`RACS` datum objects.

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

        # TODO: switch default error to values in the settings file.
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

