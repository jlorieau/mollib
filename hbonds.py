"""
MolLib functions for calculating hydrogen bonds and hydrogen positions.

   @Author:             Justin L Lorieau <jlorieau>
   @Date:               2016-08-03T12:01:01-05:00
   @Last modified by:   jlorieau
   @Last modified time: 2016-08-04T11:24:24-05:00
   @License:            Copyright 2016
"""
from mollib import Molecule
from pprint import pprint
from math import sqrt, pi, acos
import numpy as np


# This is the cutoff distance between N and O atoms to be considered a
# hydrogen bond
hydrogen_bond_cutoff = 2.5  # Angstroms

# This is the optimum NH bond length in proteins.
# 1. L. Yao, B. Vogeli, J. Ying, A. Bax, J. Am. Chem. Soc.
# 130, 16518-20 (2008).
nh_optimal = 1.023  # Angstroms


def vector_length(vector):
    "Returns the length (in A) of a vector"
    return sqrt(sum([i*i for i in vector]))


def calc_vector(atom_i, atom_j, normalize=True):
    "Returns the vector between atoms 'i' and 'j' with optional normalization."
    vec = atom_i.pos - atom_j.pos

    if normalize:
        length = vector_length(vec)
        return vec / length
    else:
        return vec


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
            # hydrogen_bond_cutoff then they're not hydrogen bonded
            ho_ij = calc_vector(h_j, o_i, normalize=False)
            r_ij = vector_length(ho_ij)
            if r_ij > hydrogen_bond_cutoff:
                continue

            # Calculate the C-O--HN angle. Vectors have to be normalized;
            # co_i is already normalized, but ho_ij isn't yet.
            ho_ij /= r_ij

            dot_product = np.dot(ho_ij, co_i)
            theta = acos(dot_product) * 180. / pi

            # Add it to the list of possible hydrogen bonds
            possible_hbonds.append(['', o_i, h_j, r_ij, theta])

    hbonds = classify_amide_hbonds(possible_hbonds)
    return hbonds


def classify_amide_hbonds(possible_hbonds):
    """Takes a list of possible hbonds from 'find_hbond_partners' and
    classifies their type.

    Classification definition from A. Grishaev, A. Bax, J. Am. Chem. Soc.
    126, 7281-92 (2004).
    """
    for hbond in possible_hbonds:
        classification, o_i, h_j, r_oh, theta = hbond

        # Check alpha-helix
        if all((h_j.residue.number - o_i.residue.number == 4,
                r_oh < 2.3,
                theta > (149.-30.) and theta < (149.+30.))):  # 149 +/- 30 deg
            hbond[0] = 'bb HN: alpha-helix'
            continue

        # Check 310-helix
        if all((h_j.residue.number - o_i.residue.number == 3,
                r_oh < 2.4,
                theta > (114.-30.) and theta < (114.+30.))):  # 114 +/- 30 deg
            hbond[0] = 'bb HN: 310-helix'
            continue

        # Check pi-helix. I guessed these theta angles since they're not in
        # the publication
        if all((h_j.residue.number - o_i.residue.number == 3,
                r_oh < 2.4,
                theta > (149.-30.) and theta < (149.+30.))):  # 49 +/- 30 deg
            hbond[0] = 'bb HN: pi-helix'
            continue

        # Check beta-sheet
        if all((r_oh < 2.2,
                theta > (155.-30.) and theta < (155.+30.))):  # 155 +/- 30 deg
            hbond[0] = 'bb HN: sheet'
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
