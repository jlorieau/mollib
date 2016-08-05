"""
MolLib functions for calculating hydrogen bonds and hydrogen positions.

   @Author:             Justin L Lorieau <jlorieau>
   @Date:               2016-08-03T12:01:01-05:00
   @Last modified by:   jlorieau
   @Last modified time: 2016-08-05T15:37:44-05:00
   @License:            Copyright 2016

TODO: add hydrogenation functions for HA, HB, and so on
TODO: add classification functions for hydrogen bonds
"""
import logging
from pprint import pprint
from math import sqrt, pi, acos
import numpy as np

try:
    from util import vector_length, calc_vector
    from mollib import Molecule
    import settings
except ImportError:
    from .util import vector_length, calc_vector
    from .mollib import Molecule
    from . import settings


def in_range(value, target, range):
    "Test whether the value is in the range of the target."
    return (target - range <= value <= target + range)


class HydrogenBond(object):
    "A basic class for storing hydrogen bond information."

    donor = None
    acceptor = None
    major_classification = ''
    minor_classification = ''
    distance = None
    angle = None
    # TODO make distsance and angle a dict so that multiple distances can be
    # measured

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            if not hasattr(self, k):
                continue
            setattr(self, k, v)

    def __repr__(self):
        s = "Hbond: "
        if self.major_classification:
            s += "{major}".format(major=self.major_classification)
        if self.minor_classification:
            s += " ({minor})".format(minor=self.minor_classification)
        s += ": "
        s = s.ljust(35)  # Even the classification column
        s += "don.({d}) - acc.({a}). ".format(d=self.donor,
                                              a=self.acceptor)
        s += 'R = {d:.1f} A, angle = {a:.0f} deg.'.format(d=self.distance,
                                                          a=self.angle)
        return s


def add_h(molecule, strip_h=True):
    """Add hydrogens to a molecule.

    :molecule:  The Molecule object to add a proton to.
    :strip_h:   If true, all hydrogens will be stripped from the molecule
                first.
    """
    if strip_h:
        mol.strip_atoms(element='H')

    def missing_message(atom_name, target_name):
        "Message to display when a proton couldn't be added."
        return '{} could not be added to {}.'.format(atom_name, target_name)

    for residue in mol.residues:
        # Pull out the relevent atoms
        n = residue.get('N', None)
        ca = residue.get('CA', None)
        c = residue.get('C', None)
        cb = residue.get('CB', None)
        c_prev = (residue.last_residue.get('C', None)
                  if residue.last_residue is not None else None)

        # Add amide protons HN -- except prolines and the first residue
        if residue.name != 'PRO' and residue.number > 1:
            r = add_one_sp2_h(molecule=mol, atom_name='HN', target_atom=n,
                              atom_1=ca, atom_2=c_prev,
                              bond_length=settings.bond_length['N-H'])
            if r is not True:  # Couldn't add atom
                logging.warning(missing_message('HN', n))

        # add HA protons -- except for Gly
        if residue.name != 'GLY':
            r = add_one_sp3_h(molecule=mol, atom_name='HA', target_atom=ca,
                              atom_1=n, atom_2=cb, atom_3=c,
                              bond_length=settings.bond_length['CA-HA'])
            if r is not True:  # Couldn't add atom
                logging.warning(missing_message('HA', ca))


# add_two_sp3_h
def add_one_sp2_h(molecule, atom_name, target_atom, atom_1, atom_2,
                  bond_length):
    """Calculate and add a single proton to an sp2 hybridized atom.

    :molecule:    The Molecule object to add a proton to.
    :atom_name:   The name of the new atom to create. ex: 'HN'
    :target_name: The Atom object to which the new proton will be added to.
    :atom_1:      The first Atom object bonded to the target_name atom.
    :atom_2:      The second Atom object bonded to the target_name atom.
    :bond_length: The length of the bond between the new proton and
                  target_atom.

    :RETURNS: True if atom was succesfully added, False if it wasn't.
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
    mol.add_atom(name=atom_name, pos=h, charge=0.0, element='H',
                 residue=target_atom.residue)
    return True


def add_one_sp3_h(molecule, atom_name, target_atom, atom_1, atom_2, atom_3,
                  bond_length):
    """Calculate and add a single proton to an sp3 hybridized atom.

    :molecule:    The Molecule object to add a proton to.
    :atom_name:   The name of the new atom to create. ex: 'HN'
    :target_name: The Atom object to which the new proton will be added to.
    :atom_1:      The first Atom object bonded to the target_name atom.
    :atom_2:      The second Atom object bonded to the target_name atom.
    :atom_3:      The third Atom object bonded to the target_name atom.
    :bond_length: The length of the bond between the new proton and
                  target_atom.

    [Returns]
        True if atom was succesfully added, False if it wasn't.
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
    mol.add_atom(name=atom_name, pos=h, charge=0.0, element='H',
                 residue=target_atom.residue)
    return True


def find_amide_hbond_partners(molecule):
    """Finds the backbone amide hydrogen bonds."""
    hbonds = find_hbond_partners(molecule=molecule,
                                 donor_name_1='O', donor_name_2='C',
                                 acceptor_name='HN')
    classify_amide_hbonds(hbonds)
    return hbonds


def find_aliphatic_hbond_partners(molecule):
    """Finds the backbone aliphatic hydrogen bonds."""
    hbonds = find_hbond_partners(molecule=molecule,
                                 donor_name_1='O', donor_name_2='C',
                                 acceptor_name='HA')
    return hbonds


def find_hbond_partners(molecule, donor_name_1, donor_name_2, acceptor_name):
    """Finds the hydrogen bond partners between donor atoms (C, O) and an
    acceptor (HN) based on distance.

    [Required]
    :molecule:      The MolLib molecule object. This function expects the
                    amide protons ('HN') to be correctly placed in
                    the molecule.
    :donor_name_1:  The first donor atom object for the hbond. ex: 'O'
    :donor_name_2:  The second donor atom object for the hbond. ex: 'C'
    :acceptor_name: The first acceptor atom object for the hbond. ex: 'HN'

    [Returns]
    :possible_hbonds:  An annotated list of possible hydrogen bonds.
    """
    # Iterate over carbonyls oxygens
    possible_hbonds = []
    for res_i in molecule.residues:
        donor_1 = res_i.get(donor_name_1, None)  # typically O
        donor_2 = res_i.get(donor_name_2, None)  # typically C
        if donor_1 is None or donor_2 is None:   # Skip if carbonyl not found
            continue

        # Construct the normalized vector for the c-o bond
        vec_i = calc_vector(donor_2, donor_1)

        # Find the nearest amide proton
        for res_j in molecule.residues:
            # Skip if acceptor (hydrogen) not found or if residue 'i' and 'j'
            # are within 1 residue from each other.
            acceptor_j = res_j.get(acceptor_name, None)
            if (acceptor_j is None or
               abs(donor_1.residue.number - acceptor_j.residue.number) < 2):
                continue

            # Calculate the donor-acceptor (HN--O) distance. If these
            # aren't within hbond_cutoff then they're not hydrogen bonded
            vec_ij = calc_vector(acceptor_j, donor_1, normalize=False)
            r_ij = vector_length(vec_ij)
            if r_ij > settings.hbond_cutoff:
                continue

            # Calculate the donor-acceptor (C-O--HN) angle. Vectors have
            # to be normalized; vec_i is already normalized, but vec_ij
            # isn't yet.
            vec_ij /= r_ij

            dot_product = np.dot(vec_ij, vec_i)
            theta = acos(dot_product) * 180. / pi

            # Add it to the list of possible hydrogen bonds
            donor = None
            acceptor = None
            major_classification = ''
            minor_classification = ''
            distance = None
            angle = None
            hbond = HydrogenBond(donor=donor_1, acceptor=acceptor_j,
                                 distance=r_ij, angle=theta)
            possible_hbonds.append(hbond)
    return possible_hbonds


def classify_amide_beta_turns(hbond):
    """Given an hbond, this function will classify beta turns.

    :hbond:     The HydrogenBond object to classify.
    """
    # Beta turns have a hydrogen bond between residue i and i+3
    if abs(hbond.acceptor.residue.number - hbond.donor.residue.number) != 3:
        return None

    # Beta turns are then defined by the ramachandran angles of the residues
    # between the residues i and i+3.
    res_0 = hbond.donor.residue
    res_1 = hbond.donor.residue.next_residue
    res_2 = hbond.acceptor.residue.last_residue
    res_3 = hbond.acceptor.residue

    phi_1, psi_1 = res_1.ramachandran_angles
    phi_2, psi_2 = res_2.ramachandran_angles

    if phi_1 is None or psi_1 is None or phi_2 is None or psi_2 is None:
        return None

    if (-80. <= phi_1 <= -40.):  # -60 +/- 20 deg
        # psi_1: -30 +/- 20 deg, phi_2: -90 +/- 20 deg
        if ((-50. <= psi_1 <= -10.) and (-110. <= phi_2 <= -70.)):
            hbond.minor_classification = "beta-turn type I"
        # psi_1: 120 +/- 20 deg, phi_2: 80 +/- 20 deg
        elif ((100. <= psi_1 <= 140.) and (60. <= phi_2 <= 100.)):
            hbond.minor_classification = "beta-turn type II"
    elif (40. <= phi_1 <= 80.):  # 60 +/- 20 deg
        # psi_1: 30 +/- 20 deg, phi_2: 90 +/- 20 deg
        if ((10. <= psi_1 <= 50.) and (70. <= phi_2 <= 110.)):
            hbond.minor_classification = "beta-turn type I'"
        # psi_1: -120 +/- 20 deg, phi_2: -80 +/- 20 deg
        elif ((-140. <= psi_1 <= -100.) and (-100. <= phi_2 <= -60.)):
            hbond.minor_classification = "beta-turn type II"
    return None


def classify_amide_hbonds(possible_hbonds):
    """Takes a list of possible hbonds from 'find_hbond_partners' and
    classifies their type.

    Classification definition from A. Grishaev, A. Bax, J. Am. Chem. Soc.
    126, 7281-92 (2004).
    """
    for hbond in possible_hbonds:
        donor = hbond.donor
        acceptor = hbond.acceptor
        distance = hbond.distance
        angle = hbond.angle

        if donor.name == 'O' and acceptor.name == 'HN':
            hbond.major_classification = 'bb HN'

        # Check beta-turn
        classify_amide_beta_turns(hbond)
        if hbond.minor_classification:
            continue

        # Check alpha-helix
        if all((acceptor.residue.number - donor.residue.number == 4,
                in_range(angle, settings.hbond_a_helix_angle,
                         settings.hbond_a_helix_angle_threshold))):
            hbond.minor_classification = 'alpha-helix'
            continue

        # Check 310-helix
        if all((acceptor.residue.number - donor.residue.number == 3,
                in_range(angle, settings.hbond_310_helix_angle,
                         settings.hbond_310_helix_angle_threshold))):
            hbond.minor_classification = '310-helix'
            continue

        # Check pi-helix. I guessed these theta angles since they're not in
        # the publication
        if all((acceptor.residue.number - donor.residue.number == 5,
                in_range(angle, settings.hbond_pi_helix_angle,
                         settings.hbond_pi_helix_angle_threshold))):
            hbond.minor_classification = 'pi-helix'
            continue

        # Check beta-sheet
        if in_range(angle, settings.hbond_beta_sheet_angle,
                    settings.hbond_beta_sheet_angle_threshold):
            hbond.minor_classification = 'beta sheet'
            continue

    # Clear weak isolated hydrogen bonds (i.e. r_oh > 2.4 Angstroms)
    # possible_hbonds = [hbond for hbond in possible_hbonds
    #                   if hbond[0] != '' and hbond[3] < 2.4]

    # Remove non-contiguous secondary structure elements
    # assert
    return possible_hbonds

if __name__ == "__main__":
    mol = Molecule('2N7J')
    add_h(mol)
#    mol.write_pdb('2MJB_H.pdb')
    hbonds = find_amide_hbond_partners(mol)
    pprint(hbonds)
    hbonds = find_aliphatic_hbond_partners(mol)
    pprint(hbonds)
