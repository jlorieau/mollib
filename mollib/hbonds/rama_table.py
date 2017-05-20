from math import exp

import mollib.core.settings
from mollib.utils import MDTable, FormattedStr
from .classify_secondary_structure import classify_residues
from . import settings


class RamaTable(MDTable):
    """A Markdown table for Ramachandran angles.

    Attributes
    ----------
    molecule: :obj:`mollib.core.Molecule`
        The molecule to create the Ramachandran table for.


    .. note:: The Rama table is created in the hbond module because the residue
              secondary structures are assigned based on their hydrogen bonds.
    """

    molecule = None

    def __init__(self, molecule, **kwargs):
        self.molecule = molecule

        # Classify the residue secondary structure
        classify_residues(molecule)

        # Setup the table. If detailed, include extra information
        if settings.rama_table_detailed:
            cols = ('Residue', 'Phi (deg)', 'Psi (deg)', 'Classification',
                    'Minor', 'E (kT) / Prob.')
        else:
            cols = ('Residue', 'Phi (deg)', 'Psi (deg)', 'Classification',
                    'E (kT) / Prob.')
        super(RamaTable, self).__init__(*cols)

        # Add the rows
        for residue in molecule.residues:
            # Skip heteroatom chains
            if '*' in residue.chain.id:
                continue

            phi, psi = residue.ramachandran_angles
            classification = residue.classification
            if isinstance(classification, tuple) and len(classification) >= 2:
                major_classification = classification[0]
                minor_classification = classification[1]
            else:
                major_classification = ''
                minor_classification = ''

            energy = getattr(residue, 'energy_ramachandran', '-')

            # If the energy has a value (float) and it is above the energy
            # cutoff, add its value to the table. Otherwise, just print
            # a '-' character, if it is within acceptable ranges.
            E_prob = ""
            if isinstance(energy, float):
                prob = exp(-1. * energy) * 100.
                E_prob = "{:>2.1f} / {:>4.1f}%".format(energy, prob)

                if energy < mollib.core.settings.energy_cutoff_good:
                    E_prob = FormattedStr(E_prob, 'green')
                elif energy < mollib.core.settings.energy_cutoff_warning:
                    E_prob = FormattedStr(E_prob, 'yellow')
                else:
                    E_prob = FormattedStr(E_prob, 'red')

            if settings.rama_table_detailed:
                self.add_row('{}.{}'.format(residue.chain.id, residue),
                             "{:>6.1f}".format(phi or 0.),
                             "{:>6.1f}".format(psi or 0.),
                             major_classification,
                             minor_classification,
                             E_prob)
            else:
                self.add_row('{}.{}'.format(residue.chain.id, residue),
                             "{:>6.1f}".format(phi or 0.),
                             "{:>6.1f}".format(psi or 0.),
                             major_classification,
                             E_prob)

