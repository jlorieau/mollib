"""
Tools to process the dipolar and chemical shift tensors for a molecule.
"""
from math import sqrt, pi
from collections import namedtuple
import logging

import numpy as np

from mollib import Molecule
from mollib.utils.interactions import interaction_label, interaction_atoms
from mollib.utils.tensors import get_Haeberlen
from mollib.utils.rotations import R
from . import settings


CSAcomponents = namedtuple('CSAcomponents',
                           'dzz dxx dyy delta eta alpha beta ref_atom1 '
                           'ref_atom2')


class Process(object):
    """Process molecules into dipolar and anisotropic chemical shift
    interactions for the A-matrix in the SVD analysis.

    The Process is a Chain-of-Responsibility pattern.

    Attributes
    ----------
    magnetic_interactions: list of dicts
        A list of magnetic interaction dicts, one for each molecule.

        - Format: ``[dict_1, dict_2, ...]`` where dict_1, etc are the magnetic
          interaction dicts for each molecule.
        - **key**: The interaction label, i.e. '14N-H' (str)
        - **value**: A tuple of the scaling constant and the 5x1 array of the
          interaction for the dipolar or chemical shift/

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
        magnetic_interactions attribute and returned.

        Returns
        -------
        magnetic_interactions: list of dicts
            A list of magnetic interaction dicts, one for each molecule.
        """

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
        value: (float, `numpy.array`)
            A tuple with the scaling constant and the array for the SVD of this
            dipole.
        """
        # Find the dipole type. This is a tuple of the form.
        dipole_type = (atom1.element, atom2.element)

        # Make the order of atoms match those in the
        # settings.default_predicted_rdcs by reversing the order if needed
        if dipole_type[::-1] in settings.default_predicted_rdcs:
            dipole_type = dipole_type[::-1]

        # Calculate or retrieve cached the static dipolar coupling constant
        if not hasattr(self, 'dcc'):
            self.dcc = {}
        if dipole_type not in self.dcc:
            # Get the gyromagnetic ratios for the atoms, based on their
            # elements.
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

        dcc = self.dcc[dipole_type]

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

        # Scale the array by the dipolar coupling. Either the dipolar coupling
        # constants from bond lengths and gyromagnetic ratios can be used
        # or pre-calculated values for each dipolar coupling are used. However,
        # pre-calculated values can only be used if they are available in
        # settings.default_predicted_rdcs.
        if (settings.calculate_from_bonds or
            dipole_type not in settings.default_predicted_rdcs):
            # Calculate from bonds
            # Convert from Angstroms to meters
            scale = dcc * r**-3 * 1.E30
        else:
            # Get the pre-calculated value
            scale = settings.default_predicted_rdcs[dipole_type]

        arr *= scale
        return (scale, arr)

    def process(self, labels=None, **kwargs):
        """Process the dipoles identified by the tensor_keys.

        Parameters
        ----------
        labels: list of str or None
            A list of interaction label strings. ex: 'A.14N-H'

        Returns
        -------
        magnetic_interactions: list of dicts
            A list of magnetic interaction dicts, one for each molecule.
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
                    for name1, name2 in settings.default_predicted_rdcs.keys():
                        # Get the correct atoms.
                        # 1. Check to see if the atom name is in this
                        #    residue.
                        # 2. Some of the atom names may include '-1' at the
                        #    end (i.e. 'C-1'), which refers to an atom in the
                        #    previous residue. Check to see if a previous
                        #    residue is specified (not None), the retrieve the
                        #    atom from the previous residue.
                        if name1 in residue:
                            a1 = residue[name1]
                        elif (name1.endswith('-1') and
                              residue.prev_residue is not None and
                              name1[:-2] in residue.prev_residue):
                            # Strip the last two characters from the name
                            # because they're '-1'. i.e.: 'C-1' becomes 'C'
                            a1 = residue.prev_residue[name1[:-2]]
                        else:
                            continue

                        if name2 in residue:
                            a2 = residue[name2]
                        elif (name2.endswith('-1') and
                              residue.prev_residue is not None and
                              name2[:-2] in residue.prev_residue):
                            a2 = residue.prev_residue[name2[:-2]]
                        else:
                            continue

                        # The atoms exist. Create an interaction label
                        key = ((a1.chain.id, a1.residue.number, a1.name),
                               (a2.chain.id, a2.residue.number, a2.name))
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
                    scale, arr = self.process_dipole(a1, a2)
                    if label in d:
                        d[label] = (d[label][0], d[label][1] + arr)
                    else:
                        d[label] = (scale, arr)

        return self.magnetic_interactions


class ProcessACS(Process):

    def get_static_tensor(self, atom):
        """Return the static tensor component dict for the given atom.

        Parameters
        ----------

        Returns
        -------
        """
        try:
            name = atom.name
            tensor_dict = settings.default_predicted_racs[name]
        except KeyError:
                msg = "A default CSA tensor for '{}' has not been specified"
                logging.error(msg.format(atom))
                return None

        # Make sure the settings have all of the fiestin
        if any(map(lambda x: x not in tensor_dict,
                   ('delta', 'eta', 'alpha', 'beta', 'ref_atom1',
                    'ref_atom2'))):
            return None

        # Return the CSA components
        tensor = get_Haeberlen(**tensor_dict)
        csa = CSAcomponents(dzz=tensor.dzz, dxx=tensor.dxx, dyy=tensor.dyy,
                            delta=tensor.delta, eta=tensor.eta,
                            alpha=tensor_dict['alpha'],
                            beta=tensor_dict['beta'],
                            ref_atom1=tensor_dict['ref_atom1'],
                            ref_atom2=tensor_dict['ref_atom2'])
        return csa

    def process_chemical_shift(self, atom):
        """Process the anisotropic chemical for the given atom.

        Parameters
        ----------
        atom: :obj:`mollib.Atom`
            The atom to calculate the chemical shift for.

        Returns
        -------
        value: (float, `numpy.array`)
            A tuple with the scaling constant and the array for the SVD of this
            CSA.
        """
        # Get the CSA tensor components for this atom
        csa = self.get_static_tensor(atom)
        residue = atom.residue
        if csa is None or residue is None:
            return None

        # get the reference atoms
        name = csa.ref_atom1
        if name.endswith('+1'):
            if residue.next_residue is None:
                return None
            ref_atom1 = residue.next_residue[name[:-2]]
        elif name.endswith('-1'):
            if residue.prev_residue is None:
                return None
            ref_atom1 = residue.prev_residue[name[:-2]]
        else:
            ref_atom1 = residue[name]

        name = csa.ref_atom2
        if name.endswith('+1'):
            if residue.next_residue is None:
                return None
            ref_atom2 = residue.next_residue[name[:-2]]
        elif name.endswith('-1'):
            if residue.prev_residue is None:
                return None
            ref_atom2 = residue.prev_residue[name[:-2]]
        else:
            ref_atom2 = residue[name]

        # Now calculate the tensor orientations and directional cosines
        vec1 = atom.pos - ref_atom1.pos
        length = np.linalg.norm(vec1)
        if length > 0.:
            vec1 /= length

        v = ref_atom1.pos - ref_atom2.pos
        length = np.linalg.norm(v)
        if length > 0.:
            v /= length

        vec2 = np.cross(vec1, v)
        length = np.linalg.norm(vec2)
        if length > 0.:
            vec2 /= length

        vec3 = np.cross(vec1, vec2)
        length = np.linalg.norm(vec3)
        if length > 0.:
            vec3 /= length

        # Compose the vectors and rotate them as needed
        vzz = vec2
        vyy = vec1
        vxx = vec3

        R_alpha = R(vzz, csa.alpha)
        vyy = np.dot(R_alpha, vyy)
        vxx = np.dot(R_alpha, vxx)

        R_beta = R(vyy, csa.beta)
        vzz = np.dot(R_beta, vzz)
        vxx = np.dot(R_beta, vxx)

        # Get the components of the tensor
        dzz = csa.dzz
        dyy = csa.dyy
        dxx = csa.dxx

        # Calculate the array
        arr = np.array((((2. / 3.) * dxx * (vxx[1] ** 2 - vxx[0] ** 2) +  # Cyy
                         (2. / 3.) * dyy * (vyy[1] ** 2 - vyy[0] ** 2) +
                         (2. / 3.) * dzz * (vzz[1] ** 2 - vzz[0] ** 2)),

                        ((2. / 3.) * dxx * (vxx[2] ** 2 - vxx[0] ** 2) +  # Czz
                         (2. / 3.) * dyy * (vyy[2] ** 2 - vyy[0] ** 2) +
                         (2. / 3.) * dzz * (vzz[2] ** 2 - vzz[0] ** 2)),

                        ((4. / 3.) * dxx * vxx[0] * vxx[1] +  # Cxy
                         (4. / 3.) * dyy * vyy[0] * vyy[1] +
                         (4. / 3.) * dzz * vzz[0] * vzz[1]),

                        ((4. / 3.) * dxx * vxx[0] * vxx[2] +  # Cxz
                         (4. / 3.) * dyy * vyy[0] * vyy[2] +
                         (4. / 3.) * dzz * vzz[0] * vzz[2]),

                        ((4. / 3.) * dxx * vxx[1] * vxx[2] +  # Cyz
                         (4. / 3.) * dyy * vyy[1] * vyy[2] +
                         (4. / 3.) * dzz * vzz[1] * vzz[2])))

        return (csa.delta * 1000., arr * 1000.)


    def process(self, labels=None, **kwargs):
        """Process the CSAs identified by the tensor_keys.

        Parameters
        ----------
        labels: list of str or None
            A list of interaction label strings. ex: 'A.14C'

        Returns
        -------
        magnetic_interactions: list of dicts
            A list of magnetic interaction dicts, one for each molecule.
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
                    for name in settings.default_predicted_racs.keys():
                        if name in residue:
                            a1 = residue[name]
                        else:
                            continue

                        # The atoms exist. Create an interaction label
                        key = ((a1.chain.id, a1.residue.number, a1.name),)
                        label = interaction_label(key)
                        labels.add(label)

        for label in labels:
            # Find the csa and procress
            for d, molecule in zip(self.magnetic_interactions, self.molecules):
                # The following function gets the atoms for the given
                # label and molecule
                atom_list = interaction_atoms(label, molecule)

                # Count the number of returned atoms. It should be 2. If it's
                # wrong, no further processing can be done.
                if len(atom_list) != 1:
                    continue

                atom = atom_list[0][0]
                return_value = self.process_chemical_shift(atom)

                if return_value is None:
                    continue

                scale, arr = return_value

                if label in d:
                    d[label] = (d[label][0], d[label][1] + arr)
                else:
                    d[label] = (scale, arr)

        return self.magnetic_interactions
