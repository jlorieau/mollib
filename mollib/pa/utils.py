"""
Utility functions for manipulating interaction labels and keys.
"""
from .data_types import RDC, RACS


def get_data_type(interaction_label):
    """Return the data type for the given interaction label."""
    return RDC if '-' in interaction_label else RACS
