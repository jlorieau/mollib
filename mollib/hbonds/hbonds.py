"""
Functions to detect hydrogen bonds in molecules.

Implementation
--------------
The current implementation:

- Finds acceptor and donor dipoles and finds acceptor--donor pairs that are
  within a distance range and angle range determined in the settings.
- Classifies backbone hydrogen bonds based on contiguous backbone-backbone
  amide torsion angles.
- The classification algorithm is more conservative than DSSP in assigning
  helices. In our implementation, all the torsion angles in the helix must
  have helical torsion angles.


.. note:: The `hbonds` attribute is added to the 'Structural Features' of the
          molecule (:obj:`mollib.Molecule`). the `hbonds` are a list of hbond
          objects (:obj:`mollib.hbond.hbond.HydrogenBond`)
"""

from collections import namedtuple
from math import pi, acos, atan2
import numpy as np

from . import settings
from .classify_hbonds import classify_hbonds
from mollib.core import (measure_distance, calc_vector, vector_length,
                         within_distance, cross)


class Dipole(namedtuple('Dipole', 'atom1 atom2')):
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

    def __repr__(self):
        if self.atom1.residue == self.atom2.residue:
            return self.atom1.__repr__() + '--' + self.atom2.name
        else:
            return self.atom1.__repr__() + '--' + self.atom2.__repr__()


class HydrogenBond(object):
    """A class for storing information on a hydrogen bond.

    Attributes
    ----------
    donor: :obj:`Dipole`
        The donor dipole of the hydrogen bond.
    acceptor: :obj:`Dipole`
        The acceptor dipole of the hydrogen bond.
    type_classification: str
        The hydrogen bond type. ex: 'bb-bb_amide'
    major_classification: str
        The classification (minor) of the hydrogen bond. ex: 'alpha-helix'
    minor_classification: str, optional
        The classification (minor) modifier. ex: 'N-term' or 'parallel'
    distances: dict
        A dict of the distances between atoms that define the hydrogen bond.
        
        - **key**: tuple of two :obj:`Atom` objects
        - **value**: the distance (in A) between the :obj:`Atom` objects
    angles: dict
        A dict of the angles between atoms that define the hydrogen bond.
        
        - **key**: tuple of three :obj:`Atom` objects
        - **value**: the angle (in deg) between the :obj:`Atom` objects
    """

    donor = None
    acceptor = None
    type_classification = ''
    major_classification = ''
    minor_classification = ''
    distances = None
    angles = None

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            if not hasattr(self, k):
                continue
            setattr(self, k, v)

    def short_repr(self):
        """The short representation of the HydrogenBond."""
        return "Hbond don.({d}) - acc.({a})".format(d=self.donor,
                                                    a=self.acceptor)

    def __repr__(self):
        """The long representation of the HydrogenBond."""
        s = self.short_repr()

        if self.type_classification:
            s += "{type}".format(type=self.type_classification)
        if self.minor_classification:
            s += " ({major})".format(major=self.major_classification)
        s += ": "
        s = s.ljust(35)  # Even the classification column

        # Add the distances, sorted by distance
        for k, v in sorted(self.distances.items(), key=lambda i: i[1]):
            # Convert strings like 'a2d2' to '{a2}{d2}'
            names = ''.join(('{', k[0:2], '}...{', k[2:4], '}'))
            names = names.format(d1=self.donor.atom1, d2=self.donor.atom2,
                                 a1=self.acceptor.atom1, a2=self.acceptor.atom2)

            s += '\n\tR({}) = {:.1f} A'.format(names, v)

        # Add the angles, sorted by angle name
        for k, v in sorted(self.angles.items(), key=lambda i: i[0]):
            s += '\n\tangle({}) = {:.1f} deg'.format(k, v)

        return s

    def __lt__(self, other):
        self_tuple = (self.donor.atom1.residue.number,
                      self.acceptor.atom1.residue.number,
                      self.type_classification,
                      self.major_classification)
        other_tuple = (other.donor.atom1.residue.number,
                       other.acceptor.atom1.residue.number,
                       other.type_classification,
                       other.major_classification)
        return self_tuple < other_tuple


def find_dipoles(molecule, donor1_elements=None, donor2_elements=None,
                 acceptor1_elements=None, acceptor2_elements=None):
    """Find all of the hydrogen bond donor and acceptor dipoles.

    Hydrogen bond dipoles are defined as:

    donor2--donor1 .... acceptor1--acceptor2

    Parameters
    ----------
    molecule: :obj:`mollib.Molecule`
        The molecule to find the dipoles in.
    donor1_elements: str, optional
        The string for the elements for the donor1 atom. The '|' string is
        supported. ex: 'HA|H'
    donor2_elements: str, optional
        The string for the elements for the donor2 atom.
    acceptor1_elements: str, optional
        The string for the elements for the acceptor1 atom.
    acceptor2_elements: str, optional
        The string for the elements for the acceptor2 atom.

    Returns
    -------
    list1, list2
        - list1. Donor dipole set
        - list2. Acceptor dipole set

    Examples
    --------
    >>> from mollib.core import Molecule
    >>> mol = Molecule('2KXA')
    >>> donor_list, acceptor_list = find_dipoles(mol)
    >>> msg = "There are {} donor dipoles and {} acceptor dipoles."
    >>> print(msg.format(len(donor_list), len(acceptor_list)))
    There are 31 donor dipoles and 31 acceptor dipoles.
    >>> donor_list[0]
    A.G1.H3--N
    >>> acceptor_list[0]
    A.G1.O--C
    """
    # Get the donor and acceptor elements values from settings if these aren't
    # specified in the function parameters,  and it splits the atom elements at
    # the '|' character.
    donor1_elements = (donor1_elements if donor1_elements is not None
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
    donor_list = list()
    acceptor_list = list()

    for atom in molecule.atoms:

        if atom.element in donor2_elements:
            # see if donor2 is bonded to an atom of type donor1_element to form
            # a dipole
            bonded_atoms = atom.bonded_atoms(sorted=True)
            for bonded_atom in bonded_atoms:

                # If a suitable dipole partner is found, create a donor diple
                if bonded_atom.element in donor1_elements:
                    dipole = Dipole(atom1=bonded_atom, atom2=atom)
                    donor_list.append(dipole)

        if atom.element in acceptor2_elements:
            # see if acceptor2 is bonded to an atom of type acceptor1_element
            # to form a dipole
            bonded_atoms = atom.bonded_atoms(sorted=True)
            for bonded_atom in bonded_atoms:

                # If a suitable dipole partner is found, create an acceptor
                # dipole
                if bonded_atom.element in acceptor1_elements:
                    dipole = Dipole(atom1=bonded_atom, atom2=atom)
                    acceptor_list.append(dipole)

    return donor_list, acceptor_list


def dipole_distances(donor_dipole, acceptor_dipole):
    """Measure the distances between atoms in a dipole.

    Parameters
    ----------
    donor_dipole: :obj:`Dipole`
        The hydrogen bond donor dipole
    acceptor_dipole: :obj:`Dipole`
        The hydrogen bond acceptor dipole
    Returns
    -------
    distance_dict: dict
        A dict with the labels (keys) and distances (values). The keys are
        formatted as a1 for acceptor1 or d2 for donor2.

          - **keys**: dipole label (str). ex: 'd1a1', 'd2a1', 'd1a2' and 'd2a2'
          - **values**: distance in A (float). ex: 3.86, 3.22, 2.10 and 2.80.

    Examples
    --------
    >>> from mollib import Molecule
    >>> mol = Molecule('2KXA')
    >>> d_dipole = Dipole(atom1=mol['A'][7]['H'], atom2=mol['A'][7]['N'])
    >>> a_dipole = Dipole(atom1=mol['A'][3]['O'], atom2=mol['A'][3]['C'])
    >>> sorted(dipole_distances(d_dipole, a_dipole).items())
    [('d1a1', 2.06), ('d1a2', 3.27), ('d2a1', 3.03), ('d2a2', 4.25)]
    """
    a1, a2 = acceptor_dipole.atom1, acceptor_dipole.atom2
    d1, d2 = donor_dipole.atom1, donor_dipole.atom2

    returned_dict = {}

    # Measure the donor1 -- acceptor1 distance
    d1a1 = measure_distance(d1, a1)

    # Do a check on the first distance. If this is too long, the rest of the
    # measurements won't be done
    if d1a1 <= settings.hbond_distance_cutoff.get('d1a1', (1.8, 3.0))[1]:
        returned_dict['d1a1'] = round(d1a1, 2)
        returned_dict['d1a2'] = round(measure_distance(d1, a2), 2)
        returned_dict['d2a1'] = round(measure_distance(d2, a1), 2)
        returned_dict['d2a2'] = round(measure_distance(d2, a2), 2)

    return returned_dict


def dipole_angles(donor_dipole, acceptor_dipole):
    """Measure the angles between the dipoles

    Parameters
    ----------
    donor_dipole: :obj:`Dipole`
        The hydrogen bond donor dipole
    acceptor_dipole: :obj:`Dipole`
        The hydrogen bond acceptor dipole

    Returns
    -------
    angle_dict: dict
        A dict with the labels (keys) and angles (values). The keys are
        formatted using the following labels.

          - **keys**: dipole label (str). ex: 'theta', 'phi'
          - **values**: angle in deg (float). ex: 168.6, 35.8

    Examples
    --------
    >>> from mollib import Molecule
    >>> mol = Molecule('2KXA')
    >>> d_dipole = Dipole(atom1=mol['A'][7]['H'], atom2=mol['A'][7]['N'])
    >>> a_dipole = Dipole(atom1=mol['A'][3]['O'], atom2=mol['A'][3]['C'])
    >>> sorted(dipole_angles(d_dipole, a_dipole).items())
    [('phi', 35.8), ('theta', 168.6)]
    """
    a1, a2 = acceptor_dipole.atom1, acceptor_dipole.atom2
    d1, d2 = donor_dipole.atom1, donor_dipole.atom2

    returned_dict = {}

    # get the x, y, z coordinate system for the acceptor.
    # The z-axis is defined by the a1-a2 vector
    z = calc_vector(a1.pos, a2.pos, normalize=True)

    # For the y-axis the *other* heavy atom bonded to a2 is needed. i.e. the
    # one that isn't the a1 atom.
    a2_bonded = [b for b in a2.bonded_heavy_atoms(sorted=True) if b != a1]
    if len(a2_bonded) < 1:
        return []
    a2_bonded = a2_bonded[0]  # Pull the first atom from the list

    v2 = calc_vector(a2.pos, a2_bonded.pos)
    y = cross(z, v2)
    y /= vector_length(y)

    # The x-axis is orthogonal to the y- and z-axes
    x = cross(y, z)
    x /= vector_length(x)

    # Calculate the theta, phi wrt the coordinate system of the a1-d1 vector
    a1d1 = calc_vector(a1.pos, d1.pos, normalize=True)

    theta = acos(np.dot(a1d1, z)) * 180. / pi
    phi = atan2(np.dot(a1d1, y), np.dot(a1d1, x)) * 180. / pi

    returned_dict["theta"] = round(theta, 1)
    returned_dict["phi"] = round(phi, 1)

    return returned_dict


def find_hbond_partners(molecule, donor1_elements=None, donor2_elements=None,
                        acceptor1_elements=None, acceptor2_elements=None,
                        dist_cutoff=None, ang_cutoff=None):
    """Finds the hydrogen bond partners between donor atoms and an
    acceptor based on distance.

    Parameters
    ----------
    molecule: :obj:`mollib.Molecule`
        The molecule to find hydrogen bonds in.
    donor1_elements: str, optional
        The string for the elements for the donor1 atom. The '|' string is
        supported. ex: 'HA|H'
    donor2_elements: str, optional
        The string for the elements for the donor2 atom.
    acceptor1_elements: str, optional
        The string for the elements for the acceptor1 atom.
    acceptor2_elements: str, optional
        The string for the elements for the acceptor2 atom.
    dist_cutoff: dict, optional
        The distance cutoff dict (see the settings file)
    ang_cutoff: dict, optional
        The angle range cutoff dict (see the settings file)


    .. note: For the donor1_elements, donor2_elements, acceptor1_elements and
             acceptor2_elements, see the documentation for :func:`find_dipoles`.

    Returns
    -------
    list of :obj:`HydrogenBond` objects

    Examples
    --------
    >>> from mollib import Molecule
    >>> mol = Molecule('2KXA')
    >>> hbonds = find_hbond_partners(mol)
    >>> print(len(hbonds))
    21
    >>> print(hbonds[0].short_repr())
    Hbond don.(A.G1.H1--N) - acc.(A.G20.O--C)
    """
    #: TODO: Implement a cache checking and force option like 'add_hydrogens'

    # Load default parameters
    dist_cutoff = (dist_cutoff if dist_cutoff is not None else
                   settings.hbond_distance_cutoff)
    ang_cutoff = (ang_cutoff if ang_cutoff is not None else
                  settings.hbond_angle_cutoff)

    # Get the acceptor and donor dipole lists
    donor_list, acceptor_list = find_dipoles(molecule,
                                             donor1_elements,
                                             donor2_elements,
                                             acceptor1_elements,
                                             acceptor2_elements)

    # Prepare dicts with the atom ids for the donor and acceptor dipoles
    # These are used to quickly create the filtered_donor and filtered_acceptor
    # sets.
    acceptor_ids = {}
    for acceptor in acceptor_list:
        id_1 = id(acceptor.atom1)
        id_2 = id(acceptor.atom2)

        s = acceptor_ids.setdefault(id_1, set())
        s.add(acceptor)

        s = acceptor_ids.setdefault(id_2, set())
        s.add(acceptor)

    donor_ids = {}
    for donor in donor_list:
        id_1 = id(donor.atom1)
        id_2 = id(donor.atom2)

        s = donor_ids.setdefault(id_1, set())
        s.add(donor)

        s = donor_ids.setdefault(id_2, set())
        s.add(donor)

    # returned list
    hbonds = []

    # For each donor-acceptor pair, detect the distance and angles with all of
    # the acceptors
    for donor_dip in donor_list:
        # Prefilter the acceptor_list using within_distance. This speeds up
        # searches of matching items.

        # First, find which distance cutoff between atoms in the dipoles is
        # largest
        largest_cutoff = None
        largest_cutoff_atoms = None
        for k, v in settings.hbond_distance_cutoff.items():
            cutoff_max = v[1]
            if largest_cutoff is None or cutoff_max > largest_cutoff:
                largest_cutoff = cutoff_max

                # split 'a1d1' to ('a1', 'd1')
                largest_cutoff_atoms = k[0:2], k[2:4]

        # If found, find all of the nearest neighbor atoms to the donor atom
        # and filter the acceptor_list based on these atoms.
        filtered_acceptor_set = acceptor_list
        if largest_cutoff:
            if 'd1' in largest_cutoff_atoms:
                donor_atom = donor_dip.atom1
            elif 'd2' in largest_cutoff_atoms:
                donor_atom = donor_dip.atom2
            else:
                donor_atom = None

            if donor_atom:
                nearest_atoms = within_distance(donor_atom,
                                                cutoff=largest_cutoff)

                # Find the donor dipoles that have one of the nearest_atoms
                # in it. We match based on atom ids because the __eq__ Atom
                # method is expensive. The new filtered_acceptor_set will only
                # contain acceptor dipoles that are within the distance cutoff
                # of the donor dipole atoms.
                filtered_acceptor_set = set()
                for nearest_atom in nearest_atoms:
                    id_1 = id(nearest_atom)
                    if id_1 in acceptor_ids:
                        filtered_acceptor_set |= acceptor_ids[id_1]

        # Find all of the acceptor dipoles that have the right distances and
        # angles to the donor dipoles.
        for acceptor_dip in filtered_acceptor_set:
            # Measure the distances between atoms in the two dipoles
            distance_dict = dipole_distances(donor_dip, acceptor_dip)

            # Check to see if the distances are all within range
            dist_test = []
            for dist_type, (min_dist, max_dist) in dist_cutoff.items():
                if dist_type not in distance_dict:
                    continue
                dist = distance_dict[dist_type]
                dist_test.append(min_dist <= dist <= max_dist)

            # This isn't a match if there were no distance constraints matched
            # or if not all the distance constraints were within range
            if len(dist_test) == 0 or not all(dist_test):
                continue
            # Measure the angles between atoms in the two dipoles
            angle_dict = dipole_angles(donor_dip, acceptor_dip)

            # Check to see if the angles are all within range
            ang_test = []
            for ang_type, (min_ang, max_ang) in ang_cutoff.items():
                if ang_type not in angle_dict:
                    continue
                ang = angle_dict[ang_type]
                ang_test.append(min_ang <= ang <= max_ang)

            # This isn't a match if there were no angles constraints matched
            # or if not all the angle constraints were within range
            if len(ang_test) == 0 or not all(ang_test):
                continue

            # This is a hydrogen bond! Create a hydrogen bond object
            hbond = HydrogenBond(donor=donor_dip, acceptor=acceptor_dip,
                                 distances=distance_dict, angles=angle_dict)

            hbonds.append(hbond)

    # Classify this hydrogen bond and add to the hbond list
    classify_hbonds(hbonds)

    # Sort the hbonds. Sorting is needed so that the output hbond listing
    # is reproducible and deterministic
    hbonds = sorted(hbonds)

    # Set the hbonds molecular parameter and return the hbonds
    molecule.set_parameter('Structural Features', 'hbonds', hbonds)
    return hbonds
