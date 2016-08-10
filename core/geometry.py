"""
Tools to measure geometries in molecules.
"""
# Author: Justin L Lorieau
# Copyright 2016

import numpy as np
from math import acos, pi, atan2
from .utils import calc_vector, vector_length


def measure_distance(atom_1, atom_2):
    """Measure the atom_1--atom_2 distance.

    Parameters
    ----------
    atom_1 : :obj:`Atom`
        The first atom.
    atom_2 : :obj:`Atom`
        The second atom.

    Returns
    -------
    distance : float
        The distance (in Angstroms)

    Examples
    --------
    >>> from mollib import Molecule
    >>> mol = Molecule('2KXA')
    >>> d = measure_distance(mol['A'][3]['CA'], mol['A'][3]['HA'])
    >>> print("{:.2f} A".format(d))
    1.08 A
    """
    v = calc_vector(atom_1, atom_2, normalize=False)
    length = vector_length(v)
    return length


def measure_angle(atom_1, atom_2, atom_3):
    """Measure the atom_1--atom_2--atom_3 angle.

    Parameters
    ----------
    atom_1: :obj:`Atom`
        The first atom of the angle to measure.
    atom_2: :obj:`Atom`
        The second atom of the angle to measure.
    atom_3: :obj:`Atom`
        The third atom of the angle to measure.

    Returns
    -------
    angle : float
        The angle (in degrees).

    Examples
    --------
    >>> from mollib import Molecule
    >>> mol = Molecule('2KXA')
    >>> gly4 = mol['A'][4]
    >>> angle = measure_angle(gly4['HA2'], gly4['CA'], gly4['HA3'])
    >>> print("{:.1f} deg".format(angle))
    109.3 deg
    """
    v1 = calc_vector(atom_2, atom_1, normalize=True)
    v2 = calc_vector(atom_2, atom_3, normalize=True)

    angle = acos(np.dot(v1, v2))
    return angle*180./pi


def measure_dihedral(atom_1, atom_2, atom_3, atom_4):
    """Measure the atom_1--atom_2--atom_3--atom_4 dihedral angle.

    Parameters
    ----------
    atom_1 : :obj:`Atom`
        The first atom.
    atom_2 : :obj:`Atom`
        The second atom.
    atom_3 : :obj:`Atom`
        The third atom.
    atom_4 : :obj:`Atom`
        The fourth atom.

    Returns
    -------
    angle : float
        The dihedral angle.

    Examples
    --------
    >>> from mollib import Molecule
    >>> mol = Molecule('2KXA')
    >>> F3 = mol['A'][3]
    >>> angle = measure_dihedral(F3['N'], F3['CA'], F3['CB'], F3['CG'])
    >>> print("{:.1f} deg".format(angle))
    -62.0 deg
    """
    # Calculate the normalized vectors.
    ab = calc_vector(atom_1, atom_2, normalize=True)
    bc = calc_vector(atom_2, atom_3, normalize=True)
    cd = calc_vector(atom_3, atom_4, normalize=True)

    # The dihedral is the angle between the plans a-b-c and b-c-d
    # The angle between these plans can be calculated from their
    # normals (cross products)
    n1 = np.cross(ab, bc)
    n2 = np.cross(bc, cd)

    # The angle between n1 and n2 can be calculated with acos. However
    # the followin atan2 relationship returns a number between 0 and
    # 2pi
    m1 = np.cross(n1, bc)
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    angle = atan2(y, x) * 180. / np.pi
    return angle


# TODO: def measure_rmsd(molecule1, molecule2, atoms=None)
