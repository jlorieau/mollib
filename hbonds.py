"""
MolLib functions for calculating hydrogen bonds and hydrogen positions.

   @Author:             Justin L Lorieau <jlorieau>
   @Date:               2016-08-03T12:01:01-05:00
   @Last modified by:   jlorieau
   @Last modified time: 2016-08-05T09:07:26-05:00
   @License:            Copyright 2016

TODO: add hydrogenation functions for HA, HB, and so on
TODO: add classification functions for hydrogen bonds
"""
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


# add_one_sp2_h(molecule, atom_name, target_atom, atom_1, atom_2, bond_length)
# add_one_sp3_h
# add_two_sp3_h
def add_backbone_hn(molecule):
    """Function to calculate and add the backbone amide protons (HN).

    Note that this function doesn't check to see if amide protons are already
    in the molecule.

    FIXME: Currently this implementation doesn't add amide protons the residue
           #1
    """
    for res_i in molecule.residues:
        # Get the adjacent N, C(i-1) and O atoms.
        n = res_i.get('N', None)
        ca = res_i.get('CA', None)
        c = (res_i.last_residue.get('C', None)
             if res_i.last_residue is not None else None)

        # If any of the atoms are None or residue is a proline, continue
        if n is None or ca is None or c is None or res_i.name == 'PRO':
            continue

        # Calculate the n-ca, n-c and bisector vectors
        nca = calc_vector(n, ca)
        nc = calc_vector(n, c)
        bisect = nca + nc
        length = vector_length(bisect)
        bisect /= length

        # calculate the hn position along the bisector
        nh_optimal = settings.bond_length['N-H']
        hn = bisect * nh_optimal + n.pos

        # Create the new 'HN' atom
        mol.add_atom(name='HN', pos=hn, charge=0.0, element='H', residue=res_i)


def find_amide_hbond_partners(molecule):
    """Finds the hydrogen bond partners between CO and N based on distance.

    [Required]
    :molecule:         The MolLib molecule object. This function expects the
                       amide protons ('HN') to be correctly placed in
                       the molecule.

    [Returns]
    :possible_hbonds:  An annotated list of possible hydrogen bonds.

    ['type classification', atom_oxygen_i, atom_hydrogen_j, oh_distance,
    coh_angle]
    """
    # Iterate over carbonyls oxygens
    possible_hbonds = []
    for res_i in molecule.residues:
        o_i = res_i.get('O', None)
        c_i = res_i.get('C', None)
        if o_i is None or c_i is None:  # Skip if carbonyl not found
            continue

        # Construct the normalized vector for the c-o bond
        co_i = calc_vector(c_i, o_i)

        # Find the nearest amide proton
        for res_j in molecule.residues:
            # Skip if hydrogen not found or if residue 'i' and 'j' are within
            # 1 residue from each other.
            h_j = res_j.get('HN', None)
            if h_j is None or abs(o_i.residue.number - h_j.residue.number) < 2:
                continue

            # Calculate the HN--O distance. If these aren't within
            # hbond_cutoff then they're not hydrogen bonded
            ho_ij = calc_vector(h_j, o_i, normalize=False)
            r_ij = vector_length(ho_ij)
            if r_ij > settings.hbond_cutoff:
                continue

            # Calculate the C-O--HN angle. Vectors have to be normalized;
            # co_i is already normalized, but ho_ij isn't yet.
            ho_ij /= r_ij

            dot_product = np.dot(ho_ij, co_i)
            theta = acos(dot_product) * 180. / pi

            # Add it to the list of possible hydrogen bonds
            donor = None
            acceptor = None
            major_classification = ''
            minor_classification = ''
            distance = None
            angle = None
            hbond = HydrogenBond(donor=o_i, acceptor=h_j, distance=r_ij,
                                 angle=theta)
            possible_hbonds.append(hbond)

    hbonds = classify_amide_hbonds(possible_hbonds)
    return hbonds


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
    mol = Molecule('2MJB')
    mol.strip_atoms(element='H')
    add_backbone_hn(mol)
    hbonds = find_amide_hbond_partners(mol)
    pprint(hbonds)
