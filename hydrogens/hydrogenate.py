"""
Functions to add hydrogens to atoms in a molecule.
"""
import logging
from math import sqrt, cos, sin

import numpy as np

from mollib.core import calc_vector, vector_length
from . import settings


def add_hydrogens(molecule, strip=True):
    """

    Parameters
    ----------
    molecule
    strip

    Returns
    -------

    """
    if strip:
        molecule.strip_atoms(element='H')

    for atom in molecule.atoms:
        if atom.element == 'H':
            continue
        topology = atom.topology
        number_hydrogens = len([i for i in topology if i.startswith('H')])
        number_heavy_atoms = len([i for i in topology if not i.startswith('H')])

        if atom.element == 'O' and number_hydrogens == 1:
            # TODO this does not work for sp3 oxygens
            add_one_sp2_h(atom, settings.bond_length['O-H'])
        elif atom.element == 'N':
            if number_hydrogens == 1:
                add_one_sp2_h(atom, settings.bond_length['N-H'])
            if number_hydrogens == 2:
                add_two_sp2_h(atom, settings.bond_length['N-H2'])
        elif atom.element == 'C':
            if number_heavy_atoms <= 2 and number_hydrogens == 1:
                add_one_sp2_h(atom, settings.bond_length['C-H'])
            if number_heavy_atoms == 3 and number_hydrogens == 1:
                add_one_sp3_h(atom, settings.bond_length['C-H'])


def add_one_sp2_h(atom, bond_length):
    """Add a single hydrogens to an sp2 hybridized atom.

    This function is useful for adding protons to add protons with a
    120-degree geometry, like backbone HNs.

    Parameters
    ----------
    atom : :obj:`atom`
        The atom to add hydrogens to.
    bond_length: float
        The length of the atom-h bond (in Angstroms).

    Returns
    -------
    bool:
        True if atom was successfully added, False if it wasn't.


    .. note:: Protons added to double-bonded atoms will respect the (E), (Z)
              assignment convention.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_angle
    >>> mol = Molecule('2KXA')
    >>> mol.strip_atoms(element='H')
    >>> F3 = mol['A'][3]
    >>> add_one_sp2_h(F3['N'], 1.0)
    True
    >>> hn = F3['H']
    >>> c_prev, n, ca = mol['A'][2]['C'], F3['N'], F3['CA']
    >>> angle1 = measure_angle(c_prev, n, hn)
    >>> angle2 = measure_angle(ca, n, hn)
    >>> print("{:.1f} {:.1f} degs".format(angle1, angle2))
    119.4 119.4 degs
    """
    bonded_heavy_atoms = atom.bonded_heavy_atoms(sorted=True)
    residue = atom.residue
    molecule = atom.molecule

    # Get the name of the hydrogen to add
    h_name = 'H' + atom.name[1:]

    # The current implementation determines the geometry based on two bonded
    # atoms.
    if len(bonded_heavy_atoms) == 2:
        # Calculate the v1, v2 and bisector vectors
        v1 = calc_vector(atom, bonded_heavy_atoms[0], normalize=True)
        v2 = calc_vector(atom, bonded_heavy_atoms[1], normalize=True)
        bisect = v1 + v2
        length = vector_length(bisect)
        bisect /= length

        normal_outplane = np.cross(v1,v2)
        normal_outplane /= vector_length(normal_outplane)

        normal_inplane = np.cross(normal_outplane, bisect)
        normal_inplane /= vector_length(normal_inplane)

        # calculate the h position along the bisector
        # cos_inplane = np.cos(in_plane * np.pi / 180.)
        # sin_inplane = np.sin(in_plane * np.pi / 180.)
        # cos_outplane = np.cos(out_plane * np.pi / 180.)
        # sin_outplane = np.sin(out_plane * np.pi / 180.)
        #
        # h = (cos_inplane * sin_outplane * normal_inplane +
        #      sin_inplane * sin_outplane * normal_outplane +
        #      cos_outplane * bisect)
        h = bisect * bond_length
        h += atom.pos

        # Create the new hydrogen atom
        molecule.add_atom(name=h_name, pos=h, element='H',
                          residue=residue)
        return True

    # If only one bonded_heavy_atom is available (like in a CO oxygen), the
    # position is inferred by looking at the the heavy atoms bonded to the
    # bonded atoms.
    elif len(bonded_heavy_atoms) == 1:
        bonded = bonded_heavy_atoms[0]

        # Find the atoms bonded to the bonded heavy atom that aren't 'atom'
        bonded_to_bonded = bonded.bonded_heavy_atoms(sorted=True)
        bonded_to_bonded = [a for a in bonded_to_bonded if a != atom]
        if len(bonded_to_bonded) < 1:
            return False
        bonded_to_bonded = bonded_to_bonded[0]

        # This algorithm will place the proton trans from the bonded_to_bonded
        # atom.
        sin60 = sqrt(3.)/2.
        cos60 = 0.5

        v1 = calc_vector(atom, bonded_heavy_atoms[0], normalize=True)
        v2 = calc_vector(bonded_heavy_atoms[0], bonded_to_bonded,
                         normalize=True)

        n_v1v2 = np.cross(v1, v2)
        n_v1v2 /= vector_length(n_v1v2)

        n = np.cross(v1, n_v1v2)
        n /= vector_length(n)

        h = (v1 * cos60 - n * sin60) * bond_length
        h += atom.pos

        # Create the new hydrogen atom
        molecule.add_atom(name=h_name, pos=h, element='H',
                          residue=residue)
        return True
    else:
        logging.warning("Number of bonded atoms "
                        "for {} is {}.".format(atom, len(bonded_heavy_atoms)))
        return False


def add_two_sp2_h(atom, bond_length):
    """Add two hydrogens to an sp2 hybridized atom.

    Parameters
    ----------
    atom : :obj:`atom`
        The atom object to add a hydrogen to.
    bond_length: float
        The length of the bond (in Angstroms) between the new proton an
        target_atom

    Returns
    -------
    bool
        True if atom was successfully added, False if it wasn't.


    .. note:: Protons added to double-bonded atoms will respect the (E), (Z)
              assignment convention.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_angle, measure_distance
    >>> mol = Molecule('2MJB')
    >>> mol.strip_atoms(element='H')
    >>> Q2 = mol['A'][2]
    >>> ne2 = Q2['NE2']
    >>> add_two_sp2_h(ne2, 1.0)
    True
    >>> h1, h2, cd, oe1 = Q2['HE21'], Q2['HE22'], Q2['CD'], Q2['OE1']
    >>> angle1 = measure_angle(h1, ne2, cd)
    >>> angle2 = measure_angle(h2, ne2, cd)
    >>> print("{:.1f} {:.1f} degs".format(angle1, angle2))
    120.0 120.0 degs
    >>> d1 = measure_distance(h1, oe1)
    >>> d2 = measure_distance(h2, oe1)
    >>> print("E:{:.1f} A, Z:{:.1f} A".format(d1, d2))  # stereoassn't
    E:3.1 A, Z:2.5 A
    """
    bonded_heavy_atoms = atom.bonded_heavy_atoms(sorted=True)
    residue = atom.residue
    molecule = atom.molecule

    # Get the name of the hydrogen to add
    h_name_1 = 'H' + atom.name[1:] + '1'  # E
    h_name_2 = 'H' + atom.name[1:] + '2'  # Z

    if len(bonded_heavy_atoms) == 1:
        bonded = bonded_heavy_atoms[0]

        # Find the atoms bonded to the bonded heavy atom that aren't 'atom'
        bonded_to_bonded = bonded.bonded_heavy_atoms(sorted=True)
        bonded_to_bonded = [a for a in bonded_to_bonded if a != atom]

        if len(bonded_to_bonded) < 1:
            return False
        bonded_to_bonded = bonded_to_bonded[0]

        # Calculate the atom--bonded_atom vector (x-axis)
        x = calc_vector(atom, bonded, normalize=True)

        # Calculate the bonded_atom--bonded_to_bonded vector. The
        # bonded_to_bonded should be sorted by mass in order to preserve
        # the E/Z notation.
        vec = calc_vector(bonded, bonded_to_bonded, normalize=True)

        # Calculate the z-axis as the orthogonal vector to these two
        z = np.cross(x, vec)
        z /= vector_length(z)

        # Calculate the y-axis as the orthogonal to the x- and y- axes.
        # The y-axis is pointing toward the bonded_to_bonded atom for
        # the 'Z' assignment
        y = np.cross(z, x)
        y /= vector_length(y)

        # The new protons are along the x- and y-axes, tilted by 60-degrees
        cos_60 = 0.5
        sin_60 = sqrt(3.) / 2.
        # The convention for E/Z assignment of Arg is the reverse for
        # Asn and Gln. Good grief. See: Markley et al., JBNMR 12, 1-23 (1998)
        if not atom.residue.name == 'ARG':
            h1_vec = x * cos_60 + y * sin_60
            h2_vec = x * cos_60 - y * sin_60
        else:
            h1_vec = x * cos_60 - y * sin_60
            h2_vec = x * cos_60 + y * sin_60

        h1 = h1_vec * bond_length + atom.pos
        h2 = h2_vec * bond_length + atom.pos

        # Add the new protons to the target_atom
        molecule.add_atom(name=h_name_1, pos=h1,
                          element='H', residue=atom.residue)
        molecule.add_atom(name=h_name_2, pos=h2,
                          element='H', residue=atom.residue)

        return True


def add_one_sp3_h(atom, bond_length):
    """Add a single hydrogen to an sp3 hybridized atom.

    This function is useful for adding methine protons, like backbone HAs.

    Parameters
    ----------
    atom : :obj:`atom`
        The atom to add a hydrogen to.
    bond_length: float
        The length of the atom-h bond (in Angstroms).

    Returns
    -------
    bool:
        True if atom was successfully added, False if it wasn't.


    Examples
    --------
    >>> from mollib.core import Molecule, measure_angle
    >>> mol = Molecule('2KXA')
    >>> mol.strip_atoms(element='H')
    >>> F3 = mol['A'][3]
    >>> n, ca, c, cb = F3['N'], F3['CA'], F3['C'], F3['CB']
    >>> add_one_sp3_h(ca, 1.0)
    True
    >>> ha = F3['HA']
    >>> angle1 = measure_angle(n, ca, ha)
    >>> angle2 = measure_angle(c, ca, ha)
    >>> angle3 = measure_angle(cb, ca, ha)
    >>> print("{:.1f} {:.1f} {:.1f} degs".format(angle1, angle2, angle3))
    108.4 108.7 109.5 degs
    """
    bonded_heavy_atoms = atom.bonded_heavy_atoms(sorted=True)
    molecule = atom.molecule

    # Get the name of the hydrogen to add
    h_name = 'H' + atom.name[1:]

    if len(bonded_heavy_atoms) == 3:
        # Calculate the plane and the plane normal formed by
        # atom_1, atom_2 and atom_3
        v1 = calc_vector(bonded_heavy_atoms[0], atom)
        v1 /= vector_length(v1)
        v2 = calc_vector(bonded_heavy_atoms[1], atom)
        v2 /= vector_length(v2)
        v3 = calc_vector(bonded_heavy_atoms[2], atom)
        v3 /= vector_length(v3)
        h_v = v1+v2+v3
        h_v /= vector_length(h_v)

        # calculate the h position along the norm
        h = -1.*h_v * bond_length + atom.pos

        # Create the new hydrogen atom
        molecule.add_atom(name=h_name, pos=h, element='H',
                          residue=atom.residue)
        return True
    else:
        msg = 'add_one_sp3_h() requires 3 heavy atoms.'
        logging.warning(msg)

