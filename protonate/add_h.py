"""
Functions to add hydrogens to molecules.
"""
# Author: Justin L Lorieau
# Copyright 2016
# TODO: add hydrogenation functions for HA, HB, and so on

import numpy as np
from mollib import settings
from mollib.core import calc_vector, vector_length


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
        "Message to display when a proton couldn't be added."
        return '{} could not be added to {}.'.format(atom_name, target_name)

    for residue in molecule.residues:
        # Pull out the relevent atoms
        n = residue.get('N', None)
        ca = residue.get('CA', None)
        c = residue.get('C', None)
        cb = residue.get('CB', None)
        c_prev = (residue.last_residue.get('C', None)
                  if residue.last_residue is not None else None)

        # Add amide protons HN -- except prolines and the first residue
        if residue.name != 'PRO' and residue.number > 1:
            r = add_one_sp2_h(molecule=molecule, atom_name='HN', target_atom=n,
                              atom_1=ca, atom_2=c_prev,
                              bond_length=settings.bond_length['N-H'])
            if r is not True:  # Couldn't add atom
                logging.warning(missing_message('HN', n))

        # add HA protons -- except for Gly
        if residue.name != 'GLY':
            r = add_one_sp3_h(molecule=molecule, atom_name='HA',
                              target_atom=ca,
                              atom_1=n, atom_2=cb, atom_3=c,
                              bond_length=settings.bond_length['CA-HA'])
            if r is not True:  # Couldn't add atom
                logging.warning(missing_message('HA', ca))


# add_two_sp3_h
def add_one_sp2_h(molecule, atom_name, target_atom, atom_1, atom_2,
                  bond_length):
    """Calculate and add a single proton to an sp2 hybridized atom.

    Parameters
    ----------
    molecule: data-type
        The Molecule object to add a proton to.
    atom_name: str
        The name of the new atom to create. ex: 'HN'
    target_name: data-type
        The Atom object to which the new proton will be added to.
    atom_1: data-type
        The first Atom object bonded to the target_name atom.
    atom_2: data-type
        The second Atom object bonded to the target_name atom.
    bond_length: float
        The length of the bond between the new proton and target_atom.

    Returns
    -------
    bool
        True if atom was succesfully added, False if it wasn't.
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
    """Calculate and add a single proton to an sp3 hybridized atom.

    Parameters
    ----------
    molecule: data-type
        The Molecule object to add a proton to.
    atom_name: data-type
        The name of the new atom to create. ex: 'HN'
    target_name: data-type
        The Atom object to which the new proton will be added to.
    atom_1: data-type
        The first Atom object bonded to the target_name atom.
    atom_2: data-type
        The second Atom object bonded to the target_name atom.
    atom_3: data-type
        The third Atom object bonded to the target_name atom.
    bond_length: float
        The length of the bond between the new proton and target_atom.

    Returns
    -------
    bool
        True if atom was successfully added, False if it wasn't.
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
