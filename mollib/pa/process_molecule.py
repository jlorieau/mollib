"""
Tools to process the dipolar and chemical shift tensors for a molecule.
"""
from math import sqrt, pi
import logging

import numpy as np

from mollib import Molecule
from mollib.utils.interactions import interaction_label, interaction_atoms
from . import settings



class Process(object):
    """Process molecules into dipolar and anisotropic chemical shift
    interactions for the SVD analysis.

    The Process is a Chain-of-Responsibility pattern.

    Attributes
    ----------
    magnetic_interactions: list of dicts
        A list of magnetic interaction dicts, one for each molecule. Magnetic
        interaction dicts have an (key) interaction string label
        and a (value) 5x1 array for the interaction.
    molecules: list of :obj:`mollib.Molecule`
        A list of molecule objects.
    _run_automatically: bool
        If True, then the process will be run automically when the Process
        parent class :meth:`process` is invoked.

    """

    magnetic_interactions = None
    molecules = None
    _run_automatically = True
    _result = None

    def __init__(self, molecules, magnetic_interactions = None):

        # Initialize the molecules attribute
        if isinstance(molecules, list):
            self.molecules = molecules
        elif isinstance(molecules, Molecule):
            self.molecules = [molecules,]

        # initialize the magnetic_interactions attribute
        if isinstance(magnetic_interactions, list):
            self.magnetic_interactions = magnetic_interactions
        else:
            self.magnetic_interactions = [dict() for m in self.molecules]

        # Init subclasses
        self._subclass_instances = []
        for sub_cls in self.__class__.__subclasses__():
            new_instance = sub_cls(molecules, magnetic_interactions)
            new_instance._run_automatically = self._run_automatically
            self._subclass_instances.append(new_instance)


    def process(self, **kwargs):
        """Process the magnetic interactions. The results are stored in the
        magnetic_interactions attribute and returned."""

        # Process all of the subclasses and store their results
        for instance in self._subclass_instances:
            result_list = instance.process(**kwargs)
            for d, result_d in zip(self.magnetic_interactions, result_list):
                d.update(result_d)

        return self.magnetic_interactions


class ProcessDipole(Process):
    """Process dipole-dipole magnetic interactions.
    """

    def process_dipole(self, atom1, atom2):
        """Process the dipole for the two given atoms.

        Parameters
        ----------
        atom1: :obj:`mollib.Atom`
            The first atom.
        atom2: :obj:`mollib.Atom`
            The second atom.

        Returns
        -------
        array: `numpy.array`
            The array for the SVD of this dipole.
        """
        # Calculate or retrieve cached the static dipolar coupling constant
        if not hasattr(self, 'dcc'):
            self.dcc = {}
        if (atom1.element, atom2.element) not in self.dcc:
            # Get the gyromagnetic ratios for the atoms, based on their elements.
            g = settings.gamma  # set the gyromagnetic ratio

            try:
                g1 = g[atom1.element]
                g2 = g[atom2.element]
            except KeyError:
                msg = ("The gyromagnetic ratio for atom {} or {} is not "
                       "specified.")
                logging.error(msg.format(atom1, atom2))
                return None

            dcc = -1. * 1.E-7 * 1.05457E-34 * g1 * g2 / (2. * pi)

            self.dcc[(atom1.element, atom2.element)] = dcc

        dcc = self.dcc[(atom1.element, atom2.element)]

        # Now calculate the bond length and directional cosines
        x, y, z = atom2.pos - atom1.pos

        r = sqrt(x * x + y * y + z * z)
        cos_x = x / r
        cos_y = y / r
        cos_z = z / r

        # Construct the array. Definition from J Biomol NMR (2010) 47:249-258.
        arr = np.array((cos_y**2 - cos_x**2,  # Cyy
                        cos_z**2 - cos_x**2,  # Czz
                        2. * cos_x * cos_y,   # Cxy
                        2. * cos_x * cos_z,   # Cxz
                        2. * cos_y * cos_z))  # Cyz

        # scale the array by the dipolar coupling. Convert from Angstroms to
        # meters
        arr *= dcc * r**-3 * 1.E30

        return arr

    def process(self, labels=None, **kwargs):
        """Process the dipoles identified by the tensor_keys.

        Parameters
        ----------
        labels: list of str or None
            A list of interaction label strings. ex: 'A.14N-H'


        Returns
        -------
        magnetic_interactions: list of dicts
            The updated magnetic_interactions
        """

        # Convert labels to a set, if it isn't already
        if isinstance(labels, list):
            labels = set(labels)
        elif isinstance(labels, set):
            pass
        else:
            labels = set()

        # If run_automically is specified as True, then it'll automatically
        # add labels (and interactions) to calculate
        if self._run_automatically:
            for molecule in self.molecules:
                for residue in molecule.residues:
                    for a1, a2 in settings.default_predicted_rdcs:
                        if a1 in residue and a2 in residue:
                            key = ((residue.chain.id, residue.number, a1),
                                   (residue.chain.id, residue.number, a2))
                            label = interaction_label(key)
                            labels.add(label)

        for label in labels:
            # Find the dipole and procress
            for d, molecule in zip(self.magnetic_interactions, self.molecules):
                # The following function gets the atoms for the given
                # label and molecule
                atom_list = interaction_atoms(label, molecule)

                # Count the number of returned atoms. It should be 2. If it's
                # wrong, no further processing can be done.
                if len(atom_list) < 1 or any([len(i) != 2 for i in atom_list]):
                    continue

                for a1, a2 in atom_list:
                    arr = self.process_dipole(a1, a2)
                    if label in d:
                        d[label] = d[label] + arr
                    else:
                        d[label] = arr

        return self.magnetic_interactions


# class ProcessNHDipole(ProcessDipole):
#     """Process the N-H dipoles (AX spin system) in the molecule."""
#
#     _run_automatically = True
#
#
#     def process(self, **kwargs):
#         """Process the NH dipoles in the molecule(s)."""
#         keys = set()
#
#         # Formulate the keys
#         for molecule in self.molecules:
#             for residue in molecule.residues:
#                 key = ((residue.chain.id, residue.number, 'N'),
#                        (residue.chain.id, residue.number, 'H'))
#                 label = interaction_label(key)
#                 keys.add(label)
#
#         # Pass the keys to the parent process. Keys that don't correspond to
#         # real atoms will be skipped
#         ProcessDipole.process(self, tensor_keys=keys)
#
#         return self.magnetic_interactions

# class ProcessACS(Process):
#
#     _run_automatically = False
#
#     def process_chemical_shift(self, atom, ref_atom1, ref_atom2):
#         """Process the anisotropic chemical for the given atom.
#
#         Parameters
#         ----------
#         atom: :obj:`mollib.Atom`
#             The atom to calculate the chemical shift fot.
#         ref_atom: :obj:`mollib.Atom`
#             The reference atom. This is a bonded heavy atom that is needed
#             to locate the PAS in reference to atom.
#
#         Returns
#         -------
#         array: `numpy.array`
#             The array for the SVD of this dipole.
#         """
#         # # Calculate or retrieve cached the static dipolar coupling constant
#         # if not hasattr(self, 'dcc'):
#         #     self.dcc = {}
#         # if (atom1.element, atom2.element) not in self.dcc:
#         #     # Get the gyromagnetic ratios for the atoms, based on their elements.
#         #     g = settings.gamma  # set the gyromagnetic ratio
#         #
#         #     try:
#         #         g1 = g[atom1.element]
#         #         g2 = g[atom2.element]
#         #     except KeyError:
#         #         msg = "The gyromagnetic ratio for atom {} or {} is not specified."
#         #         raise PartialAlignmentException(msg.format(atom1, atom2))
#         #
#         #     dcc = -1. * 1.E-7 * 1.05457E-34 * g1 * g2
#         #
#         #     self.dcc[(atom1.element, atom2.element)] = dcc
#         #
#         # dcc = self.dcc[(atom1.element, atom2.element)]
#
#         # Now calculate the tensor orientations and directional cosines
#         vec1 = atom.pos - ref_atom1.pos
#         length = np.linalg.norm(vec1)
#         if length > 0.:
#             vec1 /= length
#
#         v = ref_atom1.pos - ref_atom2.pos
#         length = np.linalg.norm(v)
#         if length > 0.:
#             v /= length
#
#         vec2 = np.cross(vec1, v)
#         length = np.linalg.norm(vec2)
#         if length > 0.:
#             vec2 /= length
#
#         vec3 = np.cross(vec1, vec2)
#         length = np.linalg.norm(vec3)
#         if length > 0.:
#             vec3 /= length
#
#
#         vzz = vec2  # vec2
#         vyy = vec1  # vec1
#         vxx = vec3  # vec3
#
#
#         dzz = 5.8   #ppm
#         dyy = -5.8
#         dxx = 0.0
#
#         arr = np.array((((2. / 3.) * dxx * (vxx[1] ** 2 - vxx[0] ** 2) +  # Cyy
#                          (2. / 3.) * dyy * (vyy[1] ** 2 - vyy[0] ** 2) +
#                          (2. / 3.) * dzz * (vzz[1] ** 2 - vzz[0] ** 2)),
#
#                         ((2. / 3.) * dxx * (vxx[2] ** 2 - vxx[0] ** 2) +  # Czz
#                          (2. / 3.) * dyy * (vyy[2] ** 2 - vyy[0] ** 2) +
#                          (2. / 3.) * dzz * (vzz[2] ** 2 - vzz[0] ** 2)),
#
#                         ((4. / 3.) * dxx * vxx[0] * vxx[1] +  # Cxy
#                          (4. / 3.) * dyy * vyy[0] * vyy[1] +
#                          (4. / 3.) * dzz * vzz[0] * vzz[1]),
#
#                         ((4. / 3.) * dxx * vxx[0] * vxx[2] +  # Cxz
#                          (4. / 3.) * dyy * vyy[0] * vyy[2] +
#                          (4. / 3.) * dzz * vzz[0] * vzz[2]),
#
#                         ((4. / 3.) * dxx * vxx[1] * vxx[2] +  # Cyz
#                          (4. / 3.) * dyy * vyy[1] * vyy[2] +
#                          (4. / 3.) * dzz * vzz[1] * vzz[2])))
#
#         return arr * 1000.
#
# class ProcessHACS(ProcessACS):
#     """Process the N-H dipoles (AX spin system) in the molecule."""
#
#     _run_automatically = True
#
#
#     def process(self, **kwargs):
#         """Process the H ACSs in the molecule(s)."""
#
#
#         for d, molecule in zip(self.magnetic_interactions, self.molecules):
#
#             for residue in molecule.residues:
#                 if 'H' not in residue:
#                     continue
#
#                 prev_residue = residue.prev_residue
#
#                 if (prev_residue is None or
#                     'N' not in residue or
#                     'C' not in prev_residue):
#                     continue
#
#                 H = residue['H']
#                 N = residue['N']
#                 C = prev_residue['C']
#
#                 key = interaction_label(residue.chain.id, residue.number, 'H')
#                 arr = self.process_chemical_shift(H, N, C)
#                 d[key] = arr
#
#         return self.magnetic_interactions
