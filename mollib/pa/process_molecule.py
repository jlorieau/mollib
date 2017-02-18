"""
Tools to process the dipolar and chemical shift tensors for a molecule.
"""
from math import sqrt

import numpy as np

from mollib import Molecule
from . import settings


class PartialAlignmentException(Exception):
    "Exception in calculating the partial alignment."
    pass


def get_subclasses(cls):
    """Generate of all sublcasses for a class."""
    subclasses = []

    for subclass in cls.__subclasses__():
        subclasses.append(subclass)
        subclasses.extend(get_subclasses(subclass))

    return subclasses


class Process(object):
    """Process molecules into dipolar and anisotropic chemical shift
    interactions for the SVD analysis.

    Attributes
    ----------
    magnetic_interactions: list of dicts
        A list of magnetic interaction dicts, one for each molecule. Magnetic
        interaction dicts have a tensor string
        id as a key and an array for the tensor interaction as a value.
    molecules: list of :obj:`mollib.Molecule`
        A list of molecule objects.
    _run_automatically: bool
        If True, then the process will be run automically when the Process
        parent class :meth:`process` is invoked.

    """

    magnetic_interactions = None
    molecules = None
    _run_automatically = False
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


    def process(self, **kwargs):
        """Process the magnetic interactions. The results are stored in the
        magnetic_interactions attribute and returned."""
        subclasses = get_subclasses(Process)

        for subclass in subclasses:
            # Only run subclasses marked as run automatically
            if not getattr(subclass, '_run_automatically', False):
                continue

            # Create the process subclass item and process
            sub = subclass(molecules=self.molecules,
                           magnetic_interactions=self.magnetic_interactions)

            sub.process()

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
                msg = "The gyromagnetic ratio for atom {} or {} is not specified."
                raise PartialAlignmentException(msg.format(atom1, atom2))

            dcc = -1. * 1.E-7 * 1.05457E-34 * g1 * g2

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


    def process(self, tensor_keys, **kwargs):
        """Process the dipoles identified by the tensor_keys.

        Parameters
        ----------
        tensor_keys: list of tuples
            A list of tuples of the form
            [((subunit, res_number, atom_name),
               subunit, res_number, atom_name))]

        Returns
        -------
        magnetic_interactions: list of dicts
            The updated magnetic_interactions
        """

        for key in tensor_keys:
            # Check that the key match a dipolar
            try:
                ((sub1, res1, name1), (sub2, res2, name2)) = key
            except ValueError:
                continue

            # Find the dipole and procress
            for d, molecule in zip(self.magnetic_interactions, self.molecules):
                try:
                    a1 = molecule[sub1][res1][name1]
                    a2 = molecule[sub2][res2][name2]
                except KeyError:
                    continue

                arr = self.process_dipole(a1, a2)
                d[key] = arr

        return self.magnetic_interactions



class ProcessNHDipole(ProcessDipole):
    """Process the N-H dipoles (AX spin system) in the molecule."""

    _run_automatically = True


    def process(self, **kwargs):
        """Process the NH dipoles in the molecule(s)."""
        keys = set()

        # Formulate the keys
        for molecule in self.molecules:
            for residue in molecule.residues:
                key = ((residue.chain.id, residue.number, 'N'),
                       (residue.chain.id, residue.number, 'H'))
                keys.add(key)

        # Pass the keys to the parent process. Keys that don't correspond to
        # real atoms will be skipped
        ProcessDipole.process(self, tensor_keys=keys)

        return self.magnetic_interactions

class ProcessACS:
    pass