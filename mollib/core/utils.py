"""
Utility functions.
"""
from math import sqrt

import re
import itertools


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


def group_by_2(iterable):
    """Group the items of iterable into a list of tuples of 2 items

    Examples
    --------
    >>> group_by_2('ABCDE')
    [('A', 'B'), ('B', 'C'), ('C', 'D'), ('D', 'E')]
    """
    return [(iterable[i], iterable[i + 1])
            for i in range(len(iterable)) if i + 1 < len(iterable)]


def filter_atoms(*atoms, **filters):
    """Test whether the given atoms are excluded by the given filters.

    Parameters
    ----------
    atoms: list of :obj:`mollib.Atom`
        Atoms to filter
    filters: dict
        A dict with various options to filter the results.
            - different: bool, optional:
              If True (default), filter out atoms that are all the same.

            - only_intra: bool, optional.
              If True, only intraresidue measurements are displayed. This
              option is mutually exclusive with residue_delta and exclude_intra.

            - exclude_intra: bool, optional.
              If True, distances within the same residue are excluded.

            - only_intra_chain: bool, optional.
              If True, only measurements within the same chain are included.

            - exclude_intra_chain: bool, optional.
              If True, measurements within the same chain are excluded.

            - residue_delta: int, optional.
              If specified, only residues with the given residue number spacing
              are reported.

            - bonded: bool, optional.
              If True, only measurements from bonded atoms will be returned.
              The bonded is only checked for sequential atoms. For example, if
              3 atoms are passed, then the following bonds are checked:
              atom1--atom2 and atom2--atom3. A bond between atom1--atom3 is not
              checked, as this would require a cyclo molecule.

    Returns
    -------
    bool
        False if the atoms pass the criteria in filters and should not be
        filtered.
        True if the atoms fail the given filters, and they should be excluded

    Examples
    --------
    >>> from mollib.core import Molecule, filter_atoms
    >>> mol = Molecule('2MUV')
    >>> a1, a2, a3 = mol['A'][23]['N'], mol['A'][23]['CA'], mol['A'][23]['C']
    >>> filter_atoms(a1, a2, a3)
    False
    >>> filter_atoms(a1, a2, a3, exclude_intra=True)
    True
    """
    different = filters.get('different', True)

    only_intra = filters.get('only_intra', False)
    exclude_intra = filters.get('exclude_intra', False)

    only_intra_chain = filters.get('only_intra_chain', False)
    exclude_intra_chain = filters.get('exclude_intra_chain', False)

    residue_delta = filters.get('residue_delta', None)
    bonded = filters.get('bonded', False)

    # Filter atoms that are all the same, if the different filter is set
    if different and all(i == j for i,j in itertools.combinations(atoms, 2)):
        return True

    # Check that all of the atoms come from the same residue, if only_intra
    if only_intra and not all(i.residue == j.residue
                             for i,j in itertools.combinations(atoms, 2)):
        return True

    # Filter if any of the atoms belong to multiple residues and
    # exclude_intra is True
    if exclude_intra and all(i.residue == j.residue
                             for i,j in itertools.combinations(atoms, 2)):
        return True

    # Filter based on chains
    if exclude_intra_chain and all(i.chain.id == j.chain.id
                             for i,j in itertools.combinations(atoms, 2)):
        return True

    if only_intra_chain and any(i.chain.id != j.chain.id
                                for i,j in itertools.combinations(atoms, 2)):
        return True

    # Calculate the residue number differences (deltas) between residues
    # This will produce a list like [1, 2, 1]
    deltas = [j.residue.number - i.residue.number
              for i,j in itertools.combinations(atoms, 2)]

    # Check that at least one set of residues are 'residue_delta' apart, if
    # residue_delta is specified. Test that none of the deltas is equal to
    # residue_delta, then the atoms should be filtered and True returned
    if residue_delta and all(d != residue_delta for d in deltas):
        return True

    # If at least one of the residues is residue_delta apart, make sure the
    # other ones are within residue_delta apart
    if residue_delta and not all(d <= residue_delta for d in deltas):
        return True

    # Filter if the sequence of atoms are  not bonded and bonded is true.
    # Using the mollib.Atom.in_topology function because it's 30% faster than
    # mollib.Atom.bonded(). Also this is last because it is quite expensive.
    if bonded and any(not j.in_topology(i) for i, j in group_by_2(atoms)):
        return True
    return False