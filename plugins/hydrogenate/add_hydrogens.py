"""
Functions to add hydrogens to molecules.
"""
# Author: Justin L Lorieau
# Copyright 2016
# TODO: add hydrogenation functions for HA, HB, and so on
# TODO: Make all of the add_h functions have customizable angles from settings.

import logging
import numpy as np
from math import sqrt
from mollib.core import calc_vector, vector_length
from .topology import aminoacids
import settings


def add_h(molecule, strip_h=True):
    """Add hydrogens to a molecule.

    Parameters
    ----------
    molecule: :obj:`Molecule`
        The Molecule object to add a proton to.
    strip_h: bool, optional
        If true, all hydrogens will be stripped from the molecule first.


    .. note:: Methylene protons are added stereospecifically such that pro-R
              hydrogens are H2 and pro-S hydrogens are H3.
    """
    if strip_h:
        molecule.strip_atoms(element='H')

    def missing_message(atom_name, target_name):
        """Message to display when a proton couldn't be added."""
        return '{} could not be added to {}.'.format(atom_name, target_name)

    first_residue = True
    for residue in molecule.residues:
        topology_found = False

        if residue.name in aminoacids:
            topology_found = True

            # Parse all of the atom classes
            for atom_name, values in aminoacids[residue.name].items():
                target_name = values['target_name']
                other_atom_names = values['other_atom_names']
                bond_length = values['bond_length']
                function_name = values['function']

                # Load the appropriate hydrogenation function
                func = globals().get(function_name)

                # Treat the first residue amino group as special
                # FIXME: This should call a function to make an sp3 NH3 or NH2
                if first_residue and atom_name == settings.amide_atom_name:
                    first_residue = False
                    continue

                # load the relevant atoms
                target_atom = residue.get(target_name, None)
                other_atoms = []
                for other_name in other_atom_names:
                    if other_name.endswith('-1'):  # Load the previous residue
                        res = residue.last_residue
                        atom = (res.get(other_name[:-2])  # Remove the '-1'
                                if res is not None else None)
                    else:
                        atom = residue.get(other_name, None)

                    other_atoms.append(atom)

                    # Debug values
                    if __debug__:
                        logging.debug({'molecule': molecule,
                                       'residue': residue,
                                       'atom_name': atom_name,
                                       'target_atom': target_atom,
                                       'other_atoms': other_atoms,
                                       'bond_length': bond_length,
                                       'function_name': function_name,
                                       'function': func})

                r = func(molecule=molecule, atom_name=atom_name,
                         target_atom=target_atom,
                         other_atoms=other_atoms,
                         bond_length=bond_length)

                if r is False:
                    logging.warning(missing_message(atom_name, target_atom))

        if not topology_found:
            logging.warning("Topology not found for {}".format(residue))


def add_one_sp2_h(molecule, atom_name, target_atom, other_atoms,
                  bond_length):
    """Add a single hydrogens to an sp2 hybridized atom.

    This function is useful for adding protons to add protons with a
    120-degree geometry, like backbone HNs.

    Parameters
    ----------
    molecule: :obj:`Molecule`
        The Molecule object to add a proton to.
    atom_name: str
        The name of the new atom to create. ex: 'HN'
    target_atom: :obj:`Atom`
        The Atom object to which the new proton will be added to.
        These are both the atoms connected to target_atom.
    other_atoms: list
        The two other :obj:`Atom` objects relevant to positioning the new
        hydrogen.
    bond_length: float
        The length of the bond (in Angstroms) between the new proton and
        target_atom.

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
    >>> n, ca, c_prev = F3['N'], F3['CA'], F3.last_residue['C']
    >>> add_one_sp2_h(mol, 'HN', n, [c_prev, ca], 1.0)
    True
    >>> hn = F3['HN']
    >>> angle1 = measure_angle(c_prev, n, hn)
    >>> angle2 = measure_angle(ca, n, hn)
    >>> print("{:.1f} {:.1f} degs".format(angle1, angle2))
    119.4 119.4 degs
    """
    atom_1, atom_2 = other_atoms

    # If any of the atoms are None, continue
    if target_atom is None or atom_1 is None or atom_2 is None:
        return False

    # Calculate the v1, v2 and bisector vectors
    v1 = calc_vector(target_atom, atom_1)
    v2 = calc_vector(target_atom, atom_2)
    bisect = v1 + v2
    length = vector_length(bisect)
    bisect /= length

    # calculate the h position along the bisector
    h = bisect * bond_length + target_atom.pos

    # Create the new hydrogen atom
    molecule.add_atom(name=atom_name, pos=h, charge=0.0, element='H',
                      residue=target_atom.residue)
    return True


def add_two_sp2_h(molecule, atom_name, target_atom, other_atoms,
                  bond_length):
    """Add two hydrogens to an sp2 hybridized atom.

    Parameters
    ----------
    molecule : :obj:`Molecule`
        The Molecule object to add a proton to.
    atom_name : str
        The name of the new atom to create. ex: 'HE2'
    target_atom : :obj:`Atom`
        The Atom object to which the new proton will be added to. ex: 'NE2' of
        Gln.
    other_atoms : list
        The two other :obj:`Atom` objects relevant to positioning the new
        hydrogen. These are two atoms: one is connected to target_atom, and the
        other is connected to the first atom to define a plane.
    bond_length: float
        The length of the bond (in Angstroms) between the new proton an
        target_atom

    Returns
    -------
    bool
        True if atom was successfully added, False if it wasn't.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_angle, measure_distance
    >>> mol = Molecule('2PTN')
    >>> mol.strip_atoms(element='H')
    >>> Q64 = mol['A'][64]
    >>> ne2, cd, oe1 = Q64['NE2'], Q64['CD'], Q64['OE1']
    >>> add_two_sp2_h(mol, 'HE2', ne2, [cd, oe1], 1.0)
    True
    >>> h1, h2 = Q64['HE21'], Q64['HE22']
    >>> angle1 = measure_angle(h1, ne2, cd)
    >>> angle2 = measure_angle(h2, ne2, cd)
    >>> print("{:.1f} {:.1f} degs".format(angle1, angle2))
    120.0 120.0 degs
    >>> d1 = measure_distance(h1, oe1)
    >>> d2 = measure_distance(h2, oe1)
    >>> print("Z:{:.1f} A, E:{:.1f} A".format(d1, d2))  # stereoassn't
    Z:2.5 A, E:3.2 A
    """
    atom_1, atom_2 = other_atoms

    # If any of the atoms are None, continue
    if target_atom is None or atom_1 is None:
        return False

    # Calculate the x-axis (target_atom--atom_1) vector
    # Calculate the z-axis as orthogonal to the plane formed by the 3
    # heave atoms.
    # Calculate the y-axis as orthogonal to the x- and y-axes.
    x = calc_vector(target_atom, atom_1, normalize=True)
    v = calc_vector(atom_1, atom_2, normalize=True)
    z = np.cross(x, v)
    length = vector_length(z)
    z /= length

    y = np.cross(x,z)
    length = vector_length(y)
    y /= length

    # The new protons are along the x- and y-axes, tilted by 60-degrees
    cos_60 = 0.5
    sin_60 = sqrt(3.)/2.
    h1_vec = x * cos_60 + y * sin_60
    h2_vec = x * cos_60 - y * sin_60

    h1 = h1_vec * bond_length + target_atom.pos
    h2 = h2_vec * bond_length + target_atom.pos

    # Add the new protons to the target_atom
    molecule.add_atom(name=atom_name + '1', pos=h1, charge=0.0,
                      element='H', residue=target_atom.residue)
    molecule.add_atom(name=atom_name + '2', pos=h2, charge=0.0,
                      element='H', residue=target_atom.residue)

    return True


def add_one_sp3_h(molecule, atom_name, target_atom, other_atoms,
                  bond_length):
    """Add a single hydrogen to an sp3 hybridized atom.

    This function is useful for adding methine protons, like backbone HAs.

    Parameters
    ----------
    molecule: :obj:`Molecule`
        The Molecule object to add a proton to.
    atom_name: str
        The name of the new atom to create. ex: 'HN'
    target_atom: :obj:`Atom`
        The Atom object to which the new proton will be added to.
        These are the three other atoms connected to target_atom.
    other_atoms: list
        The three other :obj:`Atom` objects relevant to positioning the new
        hydrogen.
    bond_length: float
        The length of the bond (in Angstrom) between the new proton and
        target_atom.

    Returns
    -------
    bool
        True if atom was successfully added, False if it wasn't.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_angle
    >>> mol = Molecule('2KXA')
    >>> mol.strip_atoms(element='H')
    >>> F3 = mol['A'][3]
    >>> n, ca, c, cb = F3['N'], F3['CA'], F3['C'], F3['CB']
    >>> add_one_sp3_h(mol, 'HA', ca, [n, cb, c], 1.0)
    True
    >>> ha = F3['HA']
    >>> angle1 = measure_angle(n, ca, ha)
    >>> angle2 = measure_angle(c, ca, ha)
    >>> angle3 = measure_angle(cb, ca, ha)
    >>> print("{:.1f} {:.1f} {:.1f} degs".format(angle1, angle2, angle3))
    109.5 108.6 108.4 degs
    """
    atom_1, atom_2, atom_3 = other_atoms

    # If any of the atoms are None, continue
    if (target_atom is None or atom_1 is None or atom_2 is None or
       atom_3 is None):
        return False

    # Calculate the plane and the plane normal formed by
    # atom_1, atom_2 and atom_3
    v1 = calc_vector(atom_1, atom_2)
    v2 = calc_vector(atom_2, atom_3)
    norm = np.cross(v1, v2)
    length = vector_length(norm)
    norm /= length

    # calculate the h position along the norm
    h = norm * bond_length + target_atom.pos

    # Create the new hydrogen atom
    molecule.add_atom(name=atom_name, pos=h, charge=0.0, element='H',
                      residue=target_atom.residue)
    return True


def add_two_sp3_h(molecule, atom_name, target_atom, other_atoms,
                  bond_length):
    """Add two hydrogens to an sp3 hybridized atom.

    This function is useful for adding methylenes, like backbone HA2/HA3 for
    glycines.

    Parameters
    ----------
    molecule: :obj:`Molecule`
        The Molecule object to add a proton to.
    atom_name: :obj:`Atom`
        The name of the new atom to create. Note that a stereospecific number
        will be added to it (2 or 3). ex: 'HA' will be 'HA2' and 'HA3'
    target_atom: :obj:`Atom`
        The Atom object to which the new proton will be added to.
        These are the other two atoms connected to target_atom.
    other_atoms: :list
        The two other :obj:`Atom` objects relevant to positioning the new
        hydrogen.
    bond_length: float
        The length of the bond (in Angstrom) between the new proton and
        target_atom.

    Returns
    -------
    bool:
        True if atom was successfully added, False if it wasn't.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_angle
    >>> mol = Molecule('2KXA')
    >>> mol.strip_atoms(element='H')
    >>> G4 = mol['A'][4]
    >>> n, ca, c = G4['N'], G4['CA'], G4['C']
    >>> add_two_sp3_h(mol, 'HA', ca, [n, c,], 1.0)
    True
    >>> ha2, ha3 = G4['HA2'], G4['HA3']
    >>> angle1 = measure_angle(ha2, ca, ha3)
    >>> angle2 = measure_angle(n, ca, ha3)
    >>> print('{:.1f} {:.1f} degs'.format(angle1, angle2))
    109.8 113.6 degs
    """
    atom_1, atom_2 = other_atoms

    # If any of the atoms are None, continue
    if target_atom is None or atom_1 is None or atom_2 is None:
        return False

    # Calculate the bisecting vector for the atom_1--target_atom--atom_2 angle
    v1 = calc_vector(target_atom, atom_1)
    v2 = calc_vector(target_atom, atom_2)
    bisect = v1 + v2
    length = vector_length(bisect)
    bisect /= length

    # Calculate the vector orthogonal to the atom_1--target_atom--atom2 plane
    orthog = np.cross(v1, bisect)
    length = vector_length(orthog)
    orthog /= length

    # Find the unit vector for hydrogens 1 and 2.
    scale_bisect = sqrt(2.)/3
    scale_orthog = 1.

    h_v1 = bisect * scale_bisect - orthog * scale_orthog
    length = vector_length(h_v1)
    h_v1 /= length

    h_v2 = bisect * scale_bisect + orthog * scale_bisect
    length = vector_length(h_v2)
    h_v2 /= length

    # The new protons point along h_v1 and h_v2 with a length of bond_length
    # from target_atom
    h1 = h_v1 * bond_length + target_atom.pos
    h2 = h_v2 * bond_length + target_atom.pos

    # Create the new hydrogen atoms
    molecule.add_atom(name=atom_name + '2', pos=h1, charge=0.0, element='H',
                      residue=target_atom.residue)
    molecule.add_atom(name=atom_name + '3', pos=h2, charge=0.0, element='H',
                      residue=target_atom.residue)
    return True


def add_three_sp3_h(molecule, atom_name, target_atom, other_atoms,
                  bond_length):
    """Add three hydrogens to an sp3 hybridized atom.

    This function is useful for adding methyls, like HBs of alanines.

    Parameters
    ----------
    molecule: :obj:`Molecule`
        The Molecule object to add a proton to.
    atom_name: :obj:`Atom`
        The name of the new atom to create. Note that a number will be added
        to it (1, 2 or 3). ex: 'HB' will be 'HB1', 'HB2' and 'HB3'
    target_atom: :obj:`Atom`
        The Atom object to which the new proton will be added to.
    other_atoms: :list
        The other two :obj:`Atom` objects relevant to positioning the new
        hydrogen. The first is bonded to target_atom and the second is bonded
        to the first.
    bond_length: float
        The length of the bond (in Angstrom) between the new proton and
        target_atom.

    Returns
    -------
    bool:
        True if atom was successfully added, False if it wasn't.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_angle, measure_distance
    >>> mol = Molecule('2KXA')
    >>> mol.strip_atoms(element='H')
    >>> A5 = mol['A'][5]
    >>> n, ca, cb = A5['N'], A5['CA'], A5['CB']
    >>> add_three_sp3_h(mol, 'HB', cb, [ca, n], 1.0)
    True
    >>> hb1, hb2, hb3 = A5['HB1'], A5['HB2'], A5['HB3']
    >>> angle1 = measure_angle(hb1, cb, ca)
    >>> angle2 = measure_angle(hb2, cb, ca)
    >>> angle3 = measure_angle(hb3, cb, ca)
    >>> print('{:.1f} {:.1f} {:.1f} degs'.format(angle1, angle2, angle3))
    109.5 109.5 109.5 degs
    >>> d1 = measure_distance(hb1, cb)
    >>> d2 = measure_distance(hb2, cb)
    >>> d3 = measure_distance(hb3, cb)
    >>> print('{:.2f} {:.2f} {:.2f} A'.format(d1, d2, d3))
    1.00 1.00 1.00 A
    """
    atom_1, atom_2 = other_atoms

    # If any of the atoms are None, continue
    if target_atom is None or atom_1 is None:
        return False

    # Calculate the target_atom--atom_1 vector
    v1 = calc_vector(target_atom, atom_1)
    length = vector_length(v1)
    v1 /= length

    # Calculate the atom_1--atom_2 vector
    v2 = calc_vector(atom_1, atom_2)
    length = vector_length(v2)
    v2 /= length

    # Calculate the normal to v1 and v2
    norm1 = np.cross(v1, v2)
    length = vector_length(norm1)
    norm1 /= length

    # Calculate the normal to norm1 and v1. v1, norm1 and norm2 form an
    # orthonormal set of vectors
    norm2 = np.cross(norm1, v1)
    length = vector_length(norm2)
    norm2 /= length

    # All three new protons are tilted away from v1. However, the first
    # is opposite to v2 to make a trans group
    c_705 = 0.33380  # cos(70.5-deg)
    s_705 = 0.94264  # sin(70.5-deg)
    h_v1 = c_705 * v1 + s_705 * norm2
    h1 = h_v1 * bond_length + target_atom.pos

    # Hydrogens 2 and 3 and first offset by 60-degrees using the normal
    # vectors v1 and v2.
    c_60 = 0.5  # cos(60-deg)
    s_60 = sqrt(3)/2.  # sin(60-deg)
    h_v2 = c_705 * v1 + s_705 * (-1. * c_60 * norm2 + s_60 * norm1)
    h2 = h_v2 * bond_length + target_atom.pos
    h_v3 = c_705 * v1 + s_705 * (-1. * c_60 * norm2 - s_60 * norm1)
    h3 = h_v3 * bond_length + target_atom.pos

    # Create the new hydrogen atoms
    molecule.add_atom(name=atom_name + '1', pos=h1, charge=0.0, element='H',
                      residue=target_atom.residue)
    molecule.add_atom(name=atom_name + '2', pos=h2, charge=0.0, element='H',
                      residue=target_atom.residue)
    molecule.add_atom(name=atom_name + '3', pos=h3, charge=0.0, element='H',
                      residue=target_atom.residue)
    return True
