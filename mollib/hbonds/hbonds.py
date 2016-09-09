"""
Functions to detect hydrogen bonds in molecules.

Implementation
--------------

Written Molecule Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Category: hbonds
    - list of :obj:`mollib.hbond.HydrogenBond` objects
"""

import logging
from collections import namedtuple
from math import pi, acos
import numpy as np

from . import settings
from mollib.core import vector_length, calc_vector, within_distance


Dipole = namedtuple('Dipole','atom1 atom2')
"""An electric dipole namedtuple.

Attributes
----------
atom1: :obj:`mollib.Atom`
    The first atom in the dipole. This is the 'H' atom in an hbond donor or
    the 'O' atom in a hbond acceptor.
atom2: :obj:`mollib.Atom`
    The second atom in the dipole. This is the 'N' atom in an hbond donor or
    the the 'C' atom in an hbond acceptor.
"""

class HydrogenBond(object):
    """A class for storing information on a hydrogen bond.

    Attributes
    ----------
    donor: :obj:`Dipole`
        The donor dipole of the hydrogen bond.
    acceptor: :obj:`Dipole`
        The acceptor dipole of the hydrogen bond.
    major_classification: str
        The classification (major) of the hydrogen bond. ex: 'backbone amide'
        or 'backbone aliphatic'
    minor_classification: str
        The classification (minor) of the hydrogen bond. ex: 'alpha-helix'
    distances: dict
        A dict of the distances between atoms that define the hydrogen bond.
        - key: tuple of two :obj:`Atom` objects
        - value: the distance (in A) between the :obj:`Atom` objects
    angles: dict
        A dict of the angles between atoms that define the hydrogen bond.
        - key: tuple of three :obj:`Atom objects
        - value: the angle (in deg) between the :obj:`Atom objects
    """

    donor = None
    acceptor = None
    major_classification = ''
    minor_classification = ''
    distances = None
    angles = None

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
        for k, v in self.distances.items():
            s += '\n\tR({})={:.1f}A'.format(k, v)
        for k, v in self.angles.items():
            s += '\n\tangle({})={:.1f}deg'.format(k, v)
        return s


def find_dipoles(molecule, donor1_elements=None, donor2_elements=None,
                 acceptor1_elements=None, acceptor2_elements=None):
    """Find all of the hydrogen bond donor and acceptor dipoles.

    Parameters
    ----------
    molecule
    donor1_elements
    donor2_elements
    acceptor1_elements
    acceptor2_elements

    Returns
    -------
    Two sets of :obj:`Dipole` objects
        set1: donor dipole set
        set2: acceptor dipole set

    Examples
    --------
    >>> from mollib.core import Molecule
    >>> mol = Molecule('2KXA')
    >>> donor_set, acceptor_set = find_dipoles(mol)
    >>> msg = "There are {} donor dipoles and {} acceptor dipoles."
    >>> print(msg.format(len(donor_set), len(acceptor_set)))
    There are 31 donor dipoles and 31 acceptor dipoles.
    >>> sorted(donor_set)[0]
    Dipole(atom1=A.A5-N, atom2=A.A5-H)
    >>> sorted(acceptor_set)[0]
    Dipole(atom1=A.A5-C, atom2=A.A5-O)
    """
    # Get the donor and acceptor elements values from settings if these aren't
    # specified in the function parameters,  and it splits the atom elements at
    # the '|' character.
    donor1_elements  = (donor1_elements if donor1_elements is not None
                       else settings.donor1_elements)
    donor2_elements = (donor2_elements if donor2_elements is not None
                       else settings.donor2_elements)
    acceptor1_elements = (acceptor1_elements if acceptor1_elements is not None
                          else settings.acceptor1_elements)
    acceptor2_elements = (acceptor2_elements if acceptor2_elements is not None
                          else settings.acceptor2_elements)

    # Splits the atom elements at the '|' character.
    donor1_elements = [i.strip() for i in donor1_elements.split('|')]
    donor2_elements = [i.strip() for i in donor2_elements.split('|')]
    acceptor1_elements = [i.strip() for i in acceptor1_elements.split('|')]
    acceptor2_elements = [i.strip() for i in acceptor2_elements.split('|')]

    # Got through the molecule atoms and find all of the donor2 and acceptor2
    # atoms. These are used because these atoms are consistently heavy atoms.
    donor_set = set()
    acceptor_set = set()

    for atom in molecule.atoms:
        if atom.element in donor2_elements:
            # see if donor2 is bonded to an atom of type donor1_element to form
            # a dipole
            bonded_atoms = atom.bonded_atoms()
            for bonded_atom in bonded_atoms:

                # If a suitable dipole partner is found, create a donor diple
                if bonded_atom.element in donor1_elements:
                    dipole = Dipole(atom1=atom, atom2=bonded_atom)
                    donor_set.add(dipole)

        elif atom.element in acceptor2_elements:
            # see if acceptor2 is bonded to an atom of type acceptor1_element
            # to form a dipole
            bonded_atoms = atom.bonded_atoms()
            for bonded_atom in bonded_atoms:

                # If a suitable dipole partner is found, create an acceptor
                # dipole
                if bonded_atom.element in acceptor1_elements:
                    dipole = Dipole(atom1=atom, atom2=bonded_atom)
                    acceptor_set.add(dipole)
        else:
            continue

    return donor_set, acceptor_set


# def dipole_distance(dipole1, dipole2, dipole_distance_cutoff=None):
#     """
#
#     Parameters
#     ----------
#     dipole1
#     dipole2
#     dipole_distance_cutoff
#
#     Returns
#     -------
#
#     """


# def find_hbond_partners(molecule,
#                         dist_cutoff=settings.hbond_distance_cutoff,
#                         ang_cutoff=settings.hbond_angle_cutoff):
#     """Finds the hydrogen bond partners between donor atoms and an
#     acceptor based on distance.
#
#     Parameters
#     ----------
#     molecule: :obj:`mollib.Molecule`
#         The molecule to find hydrogen bonds in.
#     dist_cutoff: dict, optional
#         The distance cutoff dict (see the settings file)
#     ang_cutoff: dict, optional
#         The angle range cutoff dict (see the settings file)
#
#     Returns
#     -------
#     list of :obj:`HydrogenBond` objects
#     """
#     # Find all of the donor and acceptor elements
#     elements = dict()  # A 3-member tuple with distance, donor_element,
#                        # acceptor_element
#     for cutoff_str, cutoff_dist in dist_cutoff.items():
#         donor_element, acceptor_element = cutoff_str.split('-')
#
#         rounded_dist = round(cutoff_dist, 1)  # round the distance to 1 decimal
#         elements.add((rounded_dist, donor_element, acceptor_element))
#
#     # Group the donor and acceptor elements
#     donor_elements = {i[1] for i in elements}
#     acceptor_elements = {i[2] for i in elements}
#
#     # Find all of the donor and acceptor atoms
#     donors_dict = dict()     # A dict of distance (key) and donor atoms
#     acceptors_dict = dict()  # A dict of distance (key) and acceptor atoms
#     for atom in molecule.atoms:
#         # Add the donor atom, if found
#         if atom.element in donor_elements:
#             # Find which donor set this donor belongs to
#             for distance, donor_element, _ in elements:
#                 if atom.element in donor_element:
#                     donor_set =  donors_dict.setdefault(distance, set())
#                     donor_set.add(atom)
#                     break
#         if atom.element in acceptor_elements:
#             for distance, _, acceptor_element in elements:
#                 if atom.element in acceptor_element:
#                     acceptor_set = acceptors_dict.setdefault(distance, set())
#                     acceptor_set.add(atom)
#                     break
#
#         # At this stage, we have the donors_dict with distance (key) and
#         # atoms (value) and the acceptors_dict with distance (key) and atoms
#         # (value).
#         #
#         # The next step is to go through all the donors and find an acceptor
#         # within the distance.
#         for distance, donor_set in donors_dict.items():
#             for donor_atom in donor_set:
#                 # Find the corresponding acceptor atoms
#                 if not distance in acceptors_dict:
#                     continue
#                 acceptors_set = acceptors_dict[distance]
#
#                 # Find the distances between donor and acceptor atoms
#                 # This function returns a list of (atom objects, distance)
#                 closest_acceptors = within_distance(donor_atom,
#                                         distance_cutoff=distance,
#                                         nearest_atom_selection=acceptors_set)
#
#
#
#
#         # Add the acceptor atom, if found. Each distance cutoff has to be
#         # checked
#         elif atom.element in acceptor_elements:
#             acceptor_atoms.add(atom)
#
#     # Find all of the donor and acceptor pairs
#     for donor_atom in donor_atoms:
#
#         nearest_list = within_distance(donor_atom,
#                                        distance_cutoff=1.0,
#                                        nearest_atom_selection=acceptor_atoms)
#
#
#     return
#
#     # Iterate over carbonyls oxygens
#     possible_hbonds = []
#     for res_i in molecule.residues:
#         donor_1 = res_i.get(donor_name_1, None)  # typically O
#         donor_2 = res_i.get(donor_name_2, None)  # typically C
#         if donor_1 is None or donor_2 is None:   # Skip if carbonyl not found
#             continue
#
#         # Construct the normalized vector for the c-o bond
#         vec_i = calc_vector(donor_2, donor_1)
#
#         # Find the nearest amide proton
#         for res_j in molecule.residues:
#             # Skip if acceptor (hydrogen) not found or if residue 'i' and 'j'
#             # are within 1 residue from each other.
#             acceptor_j = res_j.get(acceptor_name, None)
#             if (acceptor_j is None or
#                abs(donor_1.residue.number - acceptor_j.residue.number) < 2):
#                 continue
#
#             # Calculate the donor-acceptor (HN--O) distance. If these
#             # aren't within hbond_cutoff then they're not hydrogen bonded
#             vec_ij = calc_vector(acceptor_j, donor_1, normalize=False)
#             r_ij = vector_length(vec_ij)
#             if r_ij > hbond_cutoff:
#                 continue
#
#             # Calculate the donor-acceptor (C-O--HN) angle. Vectors have
#             # to be normalized; vec_i is already normalized, but vec_ij
#             # isn't yet.
#             vec_ij /= r_ij
#
#             dot_product = np.dot(vec_ij, vec_i)
#             theta = acos(dot_product) * 180. / pi
#
#             # Add it to the list of possible hydrogen bonds
#             distance = {'{}-{}'.format(donor_1.name, acceptor_j.name): r_ij}
#             angle = {'{}-{}-{}'.format(donor_2.name, donor_1.name,
#                                        acceptor_j.name): theta}
#             hbond = HydrogenBond(donor=donor_1, acceptor=acceptor_j,
#                                  distance=distance, angle=angle)
#
#             possible_hbonds.append(hbond)
#     return possible_hbonds
