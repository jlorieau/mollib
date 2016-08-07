"""
Functions to add hydrogens to molecules.
"""
# Author: Justin L Lorieau
# Copyright 2016
# TODO: add hydrogenation functions for HA, HB, and so on

import logging
import numpy as np
from math import sqrt
from mollib import settings
from mollib.core import calc_vector, vector_length


# Residues with a single alpha hydrogen
alpha = ['PRO', 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE',
         'LEU', 'LYS', 'MET', 'PHE', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

# Residues with a single amide hydrogen
amide = ['GLY', 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE',
         'LEU', 'LYS', 'MET', 'PHE', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

# Residues with a methine protons (except for HAs)
methines = {'ILE': [['HB', 'CB', 'CG1', 'CA', 'CG2'],],
            'LEU': [['HG', 'CG', 'CB', 'CD1', 'CD2'],],
            'VAL': [['HB', 'CB', 'CA', 'CG1', 'CG2'],],
            'THR': [['HB', 'CB', 'CA', 'CG2', 'OG1'],],
            }

# Residues with a methylene proton
methylenes = {'PRO': [['HB', 'CB', 'CA', 'CG'],
                      ['HG', 'CG', 'CB', 'CD'],
                      ['HD', 'CD', 'CG', 'N']],
              'GLY': [['HA', 'CA', 'N', 'C'],],
              'ARG': [['HB', 'CB', 'CA', 'CG'],
                      ['HG', 'CG', 'CB', 'CD'],
                      ['HD', 'CD', 'CG', 'NE']],
              'ARG': [['HB', 'CB', 'CA', 'CG'],],
              'ASN': [['HB', 'CB', 'CA', 'CG'],],
              'ASP': [['HB', 'CB', 'CA', 'CG'],],
              'CYS': [['HB', 'CB', 'CA', 'SG'],],
              'GLN': [['HB', 'CB', 'CA', 'CG'],
                      ['HG', 'CG', 'CB', 'CD'],],
              'GLU': [['HB', 'CB', 'CA', 'CG'],
                      ['HG', 'CG', 'CB', 'CD'],],
              'HIS': [['HB', 'CB', 'CA', 'CG'],],
              'ILE': [['HG1', 'CG1', 'CB', 'CD1'],],
              'LEU': [['HB', 'CB', 'CA', 'CG'],],
              'LYS': [['HB', 'CB', 'CA', 'CG'],
                      ['HG', 'CG', 'CB', 'CD'],
                      ['HD', 'CD', 'CG', 'CE'],
                      ['HE', 'CE', 'CD', 'NZ'],],
              'MET': [['HB', 'CB', 'CA', 'CG'],
                      ['HG', 'CG', 'CB', 'SD'], ],
              'PHE': [['HB', 'CB', 'CA', 'CG'],],
              'SER': [['HB', 'CB', 'CA', 'OG'],],
              'TRP': [['HB', 'CB', 'CA', 'CG'],],
              'TYR': [['HB', 'CB', 'CA', 'CG'],],
              }


def add_h(molecule, strip_h=True):
    """Add hydrogens to a molecule.

    Parameters
    ----------
    molecule: :obj:`Molecule`
        The Molecule object to add a proton to.
    strip_h: bool, optional
        If true, all hydrogens will be stripped from the molecule first.
    """
    if strip_h:
        molecule.strip_atoms(element='H')

    def missing_message(atom_name, target_name):
        """Message to display when a proton couldn't be added."""
        return '{} could not be added to {}.'.format(atom_name, target_name)

    first_residue = True
    for residue in molecule.residues:
        # Pull out the relevent atoms
        n = residue.get('N', None)
        ca = residue.get('CA', None)
        c = residue.get('C', None)
        cb = residue.get('CB', None)
        c_prev = (residue.last_residue.get('C', None)
                  if residue.last_residue is not None else None)

        # Add amide protons HN -- except prolines and the first residue
        if residue.name in amide and not first_residue:
            atom_name = settings.amide_atom_name
            r = add_one_sp2_h(molecule=molecule, atom_name=atom_name,
                              target_atom=n,
                              atom_1=ca, atom_2=c_prev,
                              bond_length=settings.bond_length['N-H'])
            if r is not True:  # Couldn't add atom
                logging.warning(missing_message('HN', n))

        # add HA protons -- except for Gly
        if residue.name in alpha:
            r = add_one_sp3_h(molecule=molecule, atom_name='HA',
                              target_atom=ca,
                              atom_1=n, atom_2=cb, atom_3=c,
                              bond_length=settings.bond_length['CA-HA'])
            if r is not True:  # Couldn't add atom
                logging.warning(missing_message('HA', ca))

        # Add methine protons
        if residue.name in methines:
            for atom_name, target_name, a_1_name, a_2_name, a_3_name in \
                methines[residue.name]:
                target_atom = residue.get(target_name, None)
                atom_1 = residue.get(a_1_name, None)
                atom_2 = residue.get(a_2_name, None)
                atom_3 = residue.get(a_3_name, None)
                r = add_one_sp3_h(molecule=molecule, atom_name=atom_name,
                                  target_atom=target_atom,
                                  atom_1=atom_1, atom_2=atom_2, atom_3=atom_3,
                                  bond_length=settings.bond_length['C-H'])
                if r is not True:
                    logging.warning(missing_message(atom_name, target_atom))

        # Add methylene protons
        if residue.name in methylenes:
            for atom_name, target_name, a_1_name, a_2_name in \
                methylenes[residue.name]:
                target_atom = residue.get(target_name, None)
                atom_1 = residue.get(a_1_name, None)
                atom_2 = residue.get(a_2_name, None)
                r = add_two_sp3_h(molecule=molecule, atom_name=atom_name,
                                  target_atom=target_atom,
                                  atom_1=atom_1, atom_2=atom_2,
                                  bond_length=settings.bond_length['C-H'])
                if r is not True:
                    logging.warning(missing_message(atom_name, target_atom))

        first_residue = False  # No longer the first residue

def add_one_sp2_h(molecule, atom_name, target_atom, atom_1, atom_2,
                  bond_length):
    """Calculate and add a single hydrogens to an sp2 hybridized atom.

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
    atom_1: :obj:`Atom`
        The first Atom object bonded to the target_name atom.
    atom_2: :obj:`Atom`
        The second Atom object bonded to the target_name atom.
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
    >>> add_one_sp2_h(mol, 'HN', n, c_prev, ca, 1.0)
    True
    >>> hn = F3['HN']
    >>> angle1 = measure_angle(c_prev, n, hn)
    >>> angle2 = measure_angle(ca, n, hn)
    >>> print("{:.1f} {:.1f} degs".format(angle1, angle2))
    119.4 119.4 degs
    """
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


def add_one_sp3_h(molecule, atom_name, target_atom, atom_1, atom_2, atom_3,
                  bond_length):
    """Calculate and add a single hydrogen to an sp3 hybridized atom.

    This function is useful for adding methine protons, like backbone HAs.

    Parameters
    ----------
    molecule: :obj:`Molecule`
        The Molecule object to add a proton to.
    atom_name: :obj:`Atom`
        The name of the new atom to create. ex: 'HN'
    target_atom: :obj:`Atom`
        The Atom object to which the new proton will be added to.
    atom_1: :obj:`Atom`
        The first Atom object bonded to the target_name atom.
    atom_2: :obj:`Atom`
        The second Atom object bonded to the target_name atom.
    atom_3: :obj:`Atom`
        The third Atom object bonded to the target_name atom.
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
    >>> F3 = mol['A'][3]
    >>> n, ca, c, cb = F3['N'], F3['CA'], F3['C'], F3['CB']
    >>> add_one_sp3_h(mol, 'HA', ca, n, cb, c, 1.0)
    True
    >>> ha = F3['HA']
    >>> angle1 = measure_angle(n, ca, ha)
    >>> angle2 = measure_angle(c, ca, ha)
    >>> angle3 = measure_angle(cb, ca, ha)
    >>> print("{:.1f} {:.1f} {:.1f} degs".format(angle1, angle2, angle3))
    109.5 108.6 108.4 degs
    """
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


def add_two_sp3_h(molecule, atom_name, target_atom, atom_1, atom_2,
                  bond_length):
    """Calculate and add two hydrogens to an sp3 hybridized atom.

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
    atom_1: :obj:`Atom`
        The first Atom object bonded to the target_name atom.
    atom_2: :obj:`Atom`
        The second Atom object bonded to the target_name atom.
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
    >>> add_two_sp3_h(mol, 'HA', ca, n, c, 1.0)
    True
    >>> ha2, ha3 = G4['HA2'], G4['HA3']
    >>> angle1 = measure_angle(ha2, ca, ha3)
    >>> angle2 = measure_angle(n, ca, ha3)
    >>> print('{:.1f} {:.1f} degs'.format(angle1, angle2))
    109.8 113.6 degs
    """
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
