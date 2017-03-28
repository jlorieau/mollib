"""
Tools to process the dipolar and chemical shift tensors for a molecule.
"""
import logging
import re
from math import sqrt, pi
from collections import namedtuple

import numpy as np

from mollib import Molecule
from mollib.utils.interactions import interaction_label, interaction_atoms
from mollib.utils.tensors import get_Haeberlen
from mollib.utils.rotations import R
from . import settings


CSAcomponents = namedtuple('CSAcomponents',
                           'dzz dxx dyy delta eta alpha beta gamma '
                           'ref_atom1 ref_atom2 order')


# A regex to process 'xyz' or 'x-yz' strings for a csa tensor vector order.
re_order = re.compile('(\+?\-?\w)')


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
            self.molecules = [molecules, ]

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
        dipole_type = "-".join((atom1.name, atom2.name))
        dipole_type_rev = "-".join((atom2.name, atom1.name))

        # Make the order of atoms match those in the
        # settings.default_predicted_rdcs by reversing the order if needed
        if dipole_type_rev in settings.default_predicted_rdcs:
            dipole_type = dipole_type_rev

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
            self.dcc[dipole_type] = dcc

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
        # Return the scaling factor and the array. The scaling factor has to
        # be multiplied by 2 because the RDC is measured from a splitting.
        # ( J+D  - J  = D )
        return (2. * scale, arr)  # TODO: Move to SVD and scale error too

    def process(self, labels=None, **kwargs):
        """Process the dipoles identified by the interaction labels.

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

        # If run_automatically is specified as True, then it'll automatically
        # add labels (and interactions) to calculate
        if self._run_automatically:
            for molecule in self.molecules:
                for residue in molecule.residues:
                    for sublabel in settings.default_predicted_rdcs.keys():
                        label = "{}.{}".format(residue.chain.id, residue.number)
                        label = label + sublabel
                        labels.add(label)

        for label in labels:
            # Find the dipole and process
            for d, molecule in zip(self.magnetic_interactions, self.molecules):
                # The following function gets the atoms for the given
                # label and molecule
                atom_list = interaction_atoms(label, molecule)

                # Count the number of returned atoms. It should be 2. If it's
                # not 2 (for a dipole), no further processing can be done.
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
                   ('delta', 'eta', 'alpha', 'beta', 'gamma', 'ref_atom1',
                    'ref_atom2', 'order'))):
            return None

        # Return the CSA components
        tensor = get_Haeberlen(**tensor_dict)
        csa = CSAcomponents(dzz=tensor.dzz, dxx=tensor.dxx, dyy=tensor.dyy,
                            delta=tensor.delta, eta=tensor.eta,
                            alpha=tensor_dict['alpha'],
                            beta=tensor_dict['beta'],
                            gamma=tensor_dict['gamma'],
                            ref_atom1=tensor_dict['ref_atom1'],
                            ref_atom2=tensor_dict['ref_atom2'],
                            order=tensor_dict['order'])
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
        if name.endswith('+1') or name.endswith('-1'):
            # Get the residue (next or prev) for the reference atom.
            ref_residue = (residue.next_residue if name.endswith('+1')
                           else residue.prev_residue)

            # Get the name of the atom for the ref atom.
            # Remove the '+1' or '-1'.
            ref_name = name[:-2]

            # set the ref_atom1, if it exists
            if ref_residue is not None and ref_name in ref_residue:
                ref_atom1 = ref_residue[ref_name]
            else:
                return None
        elif name in residue:
            # Otherwise get the atom from this residue
            ref_atom1 = residue[name]
        else:
            return None

        name = csa.ref_atom2
        if name.endswith('+1') or name.endswith('-1'):
            # Get the residue (next or prev) for the reference atom.
            ref_residue = (residue.next_residue if name.endswith('+1')
                           else residue.prev_residue)

            # Get the name of the atom for the ref atom.
            # Remove the '+1' or '-1'.
            ref_name = name[:-2]

            # set the ref_atom1, if it exists
            if ref_residue is not None and ref_name in ref_residue:
                ref_atom2 = ref_residue[ref_name]
            else:
                return None
        elif name in residue:
            # Otherwise get the atom from this residue
            ref_atom2 = residue[name]
        else:
            return None

        # Now calculate the tensor orientations and directional cosines
        vec1 = ref_atom1.pos - atom.pos
        length = np.linalg.norm(vec1)
        if length > 0.:
            vec1 /= length

        v = ref_atom1.pos - ref_atom2.pos
        length = np.linalg.norm(v)
        if length > 0.:
            v /= length

        vec2 = np.cross(v, vec1)
        length = np.linalg.norm(vec2)
        if length > 0.:
            vec2 /= length

        vec3 = np.cross(vec1, vec2)
        length = np.linalg.norm(vec3)
        if length > 0.:
            vec3 /= length

        # Compose the vectors and rotate them as needed
        # First order the vectors correctly. The order list is a list of
        # 'x', 'y' and 'z' (or '-x', '-y' or '-z) do indicate the order of
        # vectors vec1, vec2 and vec3 and whether these have to be reflected.
        order = [i for i in re_order.split(csa.order) if i]
        vectors = [vec1, vec2, vec3]
        for count, item in enumerate(order):
            if item == 'x':
                vxx = vectors[count]
            elif item == '-x':
                vxx = -1. * vectors[count]
            elif item == 'y':
                vyy = vectors[count]
            elif item == '-y':
                vyy = -1. * vectors[count]
            elif item == 'z':
                vzz = vectors[count]
            elif item == '-z':
                vzz = -1. * vectors[count]
            else:
                msg = "Could not parse CSA tensor vector order component '{}'"
                logging.error(msg.format(item))
                return None

        # vzz = vec1  # vec2
        # vyy = vec2  # vec1
        # vxx = vec3  # vec3

        # Then rotate them as needed
        R_alpha = R(vzz, csa.alpha)
        vyy = np.dot(R_alpha, vyy)
        vxx = np.dot(R_alpha, vxx)

        R_beta = R(vyy, csa.beta)
        vzz = np.dot(R_beta, vzz)
        vxx = np.dot(R_beta, vxx)

        R_gamma = R(vxx, csa.gamma)
        vzz = np.dot(R_gamma, vzz)
        vyy = np.dot(R_gamma, vyy)

        # Get the components of the tensor
        dzz = csa.dzz / csa.delta
        dyy = csa.dyy / csa.delta
        dxx = csa.dxx / csa.delta

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

        return csa.delta * 1000., arr

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

        # If run_automatically is specified as True, then it'll automatically
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

                # Count the number of returned atoms. It should be 1 for a CSA.
                # If it's not, no further processing can be done.
                if any([len(i) != 1 for i in atom_list]):
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
