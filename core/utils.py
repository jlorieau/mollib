"""
Utility functions.
"""
from math import sqrt
import re
import tempfile
import os


def vector_length(vector):
    """Returns the length (in A) of a vector"""
    return sqrt(sum([i*i for i in vector]))


def calc_vector(atom_i, atom_j, normalize=True):
    """Returns the vector between atoms 'i' and 'j' with optional
    normalization."""
    vec = atom_i.pos - atom_j.pos

    if normalize:
        length = vector_length(vec)
        return vec / length
    else:
        return vec


re_str = re.compile(r'[a-zA-Z]')
re_float = re.compile(r'-?\d+\.\d*')
re_int = re.compile(r'-?\d+')


def convert(s):
    """Convert a string 's' into either an integer, float or string.
    Strips spaces.

    >>> value = convert('  -6.065')
    >>> print("{} {}".format(value, isinstance(value, float)))
    -6.065 True
    >>> value = convert('  3 ')
    >>> print("{} {}".format(value, isinstance(value, int)))
    3 True
    >>> value = convert(' 1232 ')
    >>> print("{} {}".format(value, isinstance(value, int)))
    1232 True
    >>> value = convert('HETATM  ')
    >>> print("{} {}".format(value, isinstance(value, str)))
    HETATM True
    >>> value = convert(' HZ3  ')
    >>> print("{} {}".format(value, isinstance(value, str)))
    HZ3 True
    """
    # If the string contains any letter, return as a string
    m = re_str.search(s)
    if m:
        return str(s).strip()

    # Try extracting float
    m = re_float.search(s)
    if m:
        return float(m.group())

    # Try extracting int
    m = re_int.search(s)
    if m:
        return int(m.group())

    # All else fails, try just returning the string
    return str(s).strip()


def clear_cache():
    """Clears the temporary cache for mollib."""
    temp_path = os.path.join(tempfile.gettempdir(), 'mollib')
    for filename in os.listdir(temp_path):
        filepath = os.path.join(temp_path, filename)
        try:
            if os.path.isfile(filepath):
                os.unlink(filepath)
        except Exception as e:
            print(e)