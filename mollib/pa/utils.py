import re

from .data_types import RDC, RACS


def get_data_type(interaction_label):
    """Return the data type for the given interaction label."""
    return RDC if '-' in interaction_label else RACS

re_label_sort = re.compile(r'[A-Z]?\.?(\d+)(.+)')
def sort_key(interaction_label):
    """Generate a sort key for the given interaction_label"""
    m = re_label_sort.match(interaction_label)
    if m:
        groups = m.groups()
        return (groups[1], int(groups[0]))
    else:
        return None