"""
MolLib functions for calculating hydrogen bonds and hydrogen positions.

   @Author:             Justin L Lorieau <jlorieau>
   @Date:               2016-08-03T12:01:01-05:00
   @Last modified by:   jlorieau
   @Last modified time: 2016-08-06T06:58:57-05:00
   @License:            Copyright 2016

TODO: add classification functions for hydrogen bonds
"""
import logging
from pprint import pprint
from math import sqrt, pi, acos
import numpy as np

from mollib import settings
from mollib.protonate import add_h
from mollib.core import Molecule, vector_length, calc_vector


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
    # TODO: Make distance and angle dicts

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            if not hasattr(self, k):
                continue
            setattr(self, k, v)

    def __repr__(self):
        s = "Hbond "
        s += "don.({d}) - acc.({a}). ".format(d=self.donor,
                                              a=self.acceptor)
        if self.major_classification:
            s += "{major}".format(major=self.major_classification)
        if self.minor_classification:
            s += " ({minor})".format(minor=self.minor_classification)
        s += ": "
        s = s.ljust(35)  # Even the classification column
        for k, v in self.distance.items():
            s += '\n\tR({})={:.1f}A'.format(k, v)
        for k, v in self.angle.items():
            s += '\n\tangle({})={:.1f}deg'.format(k, v)
        return s


def find_amide_hbond_partners(molecule):
    """Finds the backbone amide hydrogen bonds."""
    hbonds = find_hbond_partners(molecule=molecule,
                                 donor_name_1='O', donor_name_2='C',
                                 acceptor_name='HN',
                                 hbond_cutoff=settings.hbond_amide_cutoff)
    classify_amide_hbonds(hbonds)
    return hbonds


def find_aliphatic_hbond_partners(molecule):
    """Finds the backbone aliphatic hydrogen bonds."""
    hbonds = find_hbond_partners(molecule=molecule,
                                 donor_name_1='O', donor_name_2='C',
                                 acceptor_name='HA',
                                 hbond_cutoff=settings.hbond_aliphatic_cutoff)
    return hbonds


def find_hbond_partners(molecule, donor_name_1, donor_name_2, acceptor_name,
                        hbond_cutoff=settings.hbond_cutoff):
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
            if r_ij > hbond_cutoff:
                continue

            # Calculate the donor-acceptor (C-O--HN) angle. Vectors have
            # to be normalized; vec_i is already normalized, but vec_ij
            # isn't yet.
            vec_ij /= r_ij

            dot_product = np.dot(vec_ij, vec_i)
            theta = acos(dot_product) * 180. / pi

            # Add it to the list of possible hydrogen bonds
            distance = {'{}-{}'.format(donor_1.name, acceptor_j.name): r_ij}
            angle = {'{}-{}-{}'.format(donor_2.name, donor_1.name,
                                       acceptor_j.name): theta}
            hbond = HydrogenBond(donor=donor_1, acceptor=acceptor_j,
                                 distance=distance, angle=angle)

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
        distance = hbond.distance['O-HN']
        angle = hbond.angle['C-O-HN']

        if donor.name == 'O' and acceptor.name == 'HN':
            hbond.major_classification = 'bb HN'

        # Check beta-turn
        classify_amide_beta_turns(hbond)
        if hbond.minor_classification:
            continue

        # Check alpha-helix
        if all((acceptor.residue.number - donor.residue.number == 4,
                in_range(angle, settings.hbond_a_helix_coh_angle,
                         settings.hbond_a_helix_angle_threshold))):
            hbond.minor_classification = 'alpha-helix'
            continue

        # Check 310-helix
        if all((acceptor.residue.number - donor.residue.number == 3,
                in_range(angle, settings.hbond_310_helix_coh_angle,
                         settings.hbond_310_helix_angle_threshold))):
            hbond.minor_classification = '310-helix'
            continue

        # Check pi-helix. I guessed these theta angles since they're not in
        # the publication
        if all((acceptor.residue.number - donor.residue.number == 5,
                in_range(angle, settings.hbond_pi_helix_coh_angle,
                         settings.hbond_pi_helix_angle_threshold))):
            hbond.minor_classification = 'pi-helix'
            continue

        # Check beta-sheet
        if in_range(angle, settings.hbond_beta_sheet_coh_angle,
                    settings.hbond_beta_sheet_angle_threshold):
            hbond.minor_classification = 'beta sheet'
            continue

    # Clear weak isolated hydrogen bonds (i.e. r_oh > 2.4 Angstroms)
    # possible_hbonds = [hbond for hbond in possible_hbonds
    #                   if hbond[0] != '' and hbond[3] < 2.4]

    # Remove non-contiguous secondary structure elements
    # assert
    return possible_hbonds
