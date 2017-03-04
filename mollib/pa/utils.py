import re

from mollib.utils.interactions import interaction_key
from .data_types import RDC, RACS


def get_data_type(interaction_label):
    """Return the data type for the given interaction label."""
    return RDC if '-' in interaction_label else RACS


def sort_key(interaction_label):
    """Generate a sort key for the given interaction_label.

    Examples
    --------
    >>> sort_key('14N-H')
    ('N-H', 'A', 14)
    >>> sort_key('13C')
    ('C', 'A', 13)
    >>> sort_key('B.14N-13C')
    ('N-C', 'B', 14)
    >>> sort_key('B.35CA-HA2')
    ('CA-HA', 'B', 35)
    """
    # Generate a key from the label
    key = interaction_key(interaction_label)

    # The interaction type is from the atom names
    atom_names = [i[2].translate(None, '0123456789#')  # strip numbers/#
                  for i in key]
    interaction_type = '-'.join(atom_names)

    # Use the first atom's chain id residue number as the sort key.
    chain_id = key[0][0]
    res_number = key[0][1]

    return (interaction_type, chain_id, res_number)

