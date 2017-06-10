"""
Tools to process the dipolar and chemical shift tensors for a molecule.
"""
import logging
import re
from math import sqrt, pi
from collections import namedtuple

import numpy as np

from mollib import Molecule
from mollib.core import cross, vector_length
from mollib.utils.interactions import (interaction_label, interaction_atoms,
                                       interaction_type)
from mollib.utils.tensors import get_Haeberlen
from mollib.utils.rotations import R
from . import logs
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


# TODO: Add tests for calculating the scaling factor from bonds or
#  pre-calculated values.
class ProcessDipole(Process):
    """Process dipole-dipole magnetic interactions.
    """

    def __init__(self, *args, **kwargs):
        self._pre_factors = {}
        super(ProcessDipole, self).__init__(*args, **kwargs)

    def get_dcc(self, label, atom1, atom2, distance=None):
        """Return the dipolar coupling constant (dcc) between atom1 and atom2.
        
        Parameters
        ----------
        label: str
            The interaction label.
        atom1: :obj:`mollib.Atom`
            The first atom.
        atom2: :obj:`mollib.Atom`
            The second atom.
        distance: float (optional)
            If specified, use the given distance (in A) between atom1 and
            atom2 in calculating the dc.
            
        Returns
        -------
        float or None
            The dipolar coupling constant (in Hz), if it was found.
            None if a dipolar coupling constant could not be calculated.
        """
        # Try first to see if it's a methyl group
        if settings.project_methyls:
            methyl_atom = self.which_methyl(atom1, atom2)
            if methyl_atom is not None:
                return (settings.default_predicted_rdcs['CA-HA'] *
                        settings.methyl_order_parameter)

        # At this point, no methyl was found. First, determine the type of
        # dipolar coupling.
        dipole_type = interaction_type(label)

        # Next, if pre-calculated scaling constants are allowed, see if one of
        # these is available
        if (not settings.calculate_from_bonds and
           dipole_type in settings.default_predicted_rdcs):
            return settings.default_predicted_rdcs[dipole_type]

        # Now, either the dipolar coupling should be calculated from bond
        # lengths or the value was not available in the pre-calculated values.
        # First, the dipolar pre-factor from the gyromagnetic ratio should be
        # calculated.
        if dipole_type not in self._pre_factors:
            # Get the gyromagnetic ratios for the atoms, based on their
            # elements.
            g = settings.gamma  # set the gyromagnetic ratio

            g1 = g[atom1.element] if atom1.element in g else None
            g2 = g[atom2.element] if atom2.element in g else None

            # Produce an error message if one of the gyromagnetic ratios
            # could not be found
            for g, atom in ((g1, atom1), (g2, atom2)):
                if g is not None:
                    continue
                msg = "The gyromagnetic ratio for atom {} is not specified."
                msg = msg.format(atom)
                if msg not in logs.errors:
                    logging.error(msg)
                    logs.errors.add(msg)
                return None

            pre_factor = -1. * 1.E-7 * 1.05457E-34 * g1 * g2 / (2. * pi)
            self._pre_factors[dipole_type] = pre_factor

        pre_factor = self._pre_factors[dipole_type]

        # Calculate the inter-atomic distance if needed
        if distance is None:
            x, y, z = atom2.pos - atom1.pos
            distance = sqrt(x * x + y * y + z * z)

        return pre_factor * distance ** -3 * 1.E30

    def which_methyl(self, atom1, atom2):
        """Determine whether one of the atoms corresponds to the carbon of a
        methyl group. 
        
        Parameters
        ----------
        atom1: :obj:`mollib.Atom`
            The first atom.
        atom2: :obj:`mollib.Atom`
            The second atom.

        Returns
        -------
        :obj:`mollib.Atom` or None
            If one of the atoms is a carbon for a methyl group, 
            return this atom.
            If neither is a methyl group, return None.
        """
        # Count the number of 'H' atoms for each carbon to see if one is
        # a methyl group
        for atom in (atom1, atom2):
            # To be a methyl, the atom has to be a carbon.
            if atom.element != 'C' and atom.element != '13C':
                continue

            bonded = atom.bonded_atoms(sorted=False)
            h_atoms = [i for i in bonded if i.element == 'H']
            if len(h_atoms) == 3:
                return atom

        # No methyl atom was found.
        return None

    def process_dipole(self, label, atom1, atom2):
        """Process the dipole for the two given atoms.

        Parameters
        ----------
        label: str
            The interaction label.
        atom1: :obj:`mollib.Atom`
            The first atom.
        atom2: :obj:`mollib.Atom`
            The second atom.

        Returns
        -------
        value: (float, `numpy.array`) or None
            - A tuple with the scaling constant and the array for the SVD of
              this dipole.
            - None is returned if the dipole could not be calculated. This can
              happen, for example, if one or both of the gyromagnetic ratios
              could not be found 
        """
        # Calculate the bond length and directional cosines
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

        # Calculate the scale for the the array of the dipolar coupling.
        scale = self.get_dcc(label, atom1, atom2, r)
        if scale is None:
            return None

        # Return the scaling factor and the array. The scaling factor has to
        # be multiplied by 2 because the RDC is measured from a splitting.
        # ( J+D  - J  = D )
        return 2. * scale, arr  # TODO: Move to SVD and scale error too

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
        if isinstance(labels, set):
            pass
        elif hasattr(labels, '__iter__'):
            labels = set(labels)
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

                # Process each dipole in the atom_list
                for a1, a2 in atom_list:
                    return_value = self.process_dipole(label, a1, a2)
                    if return_value is None:
                        continue
                    scale, arr = return_value

                    if label in d:
                        d[label] = (d[label][0], d[label][1] + arr)
                    else:
                        d[label] = (scale, arr)

        return self.magnetic_interactions


class ProcessACS(Process):
    """Process anisotropic chemical shift magnetic interactions.
    """

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
                msg = msg.format(atom)
                if msg not in logs.errors:
                    logging.error(msg)
                    logs.errors.add(msg)
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
        length = vector_length(vec1)
        if length > 0.:
            vec1 /= length

        v = ref_atom1.pos - ref_atom2.pos
        length = vector_length(v)
        if length > 0.:
            v /= length

        vec2 = cross(v, vec1)
        length = vector_length(vec2)
        if length > 0.:
            vec2 /= length

        vec3 = cross(vec1, vec2)
        length = vector_length(vec3)
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
                msg = msg.format(item)
                if msg not in logs.errors:
                    logging.error(msg)
                    logs.errors.add(msg)
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

                # Count the number of returned atoms. There should be atoms
                # returned to process further. The atoms lists shoud also each
                # have 1 item (1 atom) for a CSA interaction.
                # If it's not, no further processing can be done.
                if len(atom_list) < 1 or any([len(i) != 1 for i in atom_list]):
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
