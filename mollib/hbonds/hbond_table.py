from math import exp

import mollib.core.settings
from mollib.utils import MDTable, FormattedStr
from . import settings


class HBondTable(MDTable):
    """A Markdown table for hydrogen bond data.

    Attributes
    ----------
    hbonds: list of :obj:`mollib.hbonds.HydrogenBond`
        The list of hydrogen bonds.
    detailed: bool, optional
        If specified in the settings or the constructor, a detailed HBond
        table will be rendered
    sort_type: bool, optional
        If specified, the hbonds will be sorted by type in the table.
    """
    hbonds = None
    detailed = None
    sort_type = None

    def __init__(self, hbonds, **kwargs):

        # Default settings
        self.detailed = settings.hbond_table_detailed
        self.sort_type = settings.hbond_table_sort_type

        # Copy over settings specified in the constructor.
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

        # If specified, sort the hbonds
        if self.sort_type:
            self.hbonds = sorted(hbonds,
                                 key=lambda hb: (hb.major_classification,
                                                 hb.minor_classification,
                                                 hb.donor.atom2.residue.number))
        else:
            self.hbonds = hbonds

        # Setup the table headers and parent table
        if self.detailed:
            cols = ('Num', 'Donor', 'Acceptor', 'Parameter', 'Value')
        else:
            cols = ('Num', 'Donor', 'Acceptor', 'Classification',
                    'E (kT) / Prob.')
        super(HBondTable, self).__init__(*cols)

        # Add rows to the table for the hbonds
        for count, hb in enumerate(self.hbonds, 1):
            donor = hb.donor
            acceptor = hb.acceptor
            classification = ('{}/{}'.format(hb.major_classification,
                                             hb.minor_classification)
                              if hb.minor_classification else
                              '{}'.format(hb.major_classification))
            classification = "{} ({})".format(hb.type_classification,
                                              classification)
            energy = getattr(hb, 'energy_hbond', '-')
            if isinstance(energy, float):
                prob = exp(-1. * energy) * 100.
                E_prob = "{:>2.1f} / {:>4.1f}%".format(energy, prob)

                if energy < mollib.core.settings.energy_cutoff_good:
                    E_prob = FormattedStr(E_prob, 'green')
                elif energy < mollib.core.settings.energy_cutoff_warning:
                    E_prob = FormattedStr(E_prob, 'yellow')
                else:
                    E_prob = FormattedStr(E_prob, 'red')
            else:
                E_prob = '-'

            if not self.detailed:
                self.add_row(count, donor, acceptor, classification, E_prob)
            else:
                self.add_row(count, donor, acceptor, 'Classification',
                             classification)
                self.add_row('', '', '', 'E (kT) / Prob.', E_prob)

                # Now add the rows for all of the distances
                dists = [(atoms, d)
                         for atoms, d in sorted(hb.distances.items(),
                                                key=lambda i: i[1])]

                # Convert the distance names and distances to strings
                # Add the atom names
                for i in range(len(dists)):
                    a = dists[i][0]
                    dist = dists[i][1]
                    name = ''.join(('{', a[0:2], '}...{', a[2:4], '}'))
                    name = name.format(a1=hb.acceptor.atom1,
                                       a2=hb.acceptor.atom2,
                                       d1=hb.donor.atom1,
                                       d2=hb.donor.atom2)

                    dist = '{:3.2f}'.format(dist)
                    self.add_row('', '', '', name, dist)

                # Get the dipole angles
                angs = [(name, a) for name, a in sorted(hb.angles.items(),
                                                        key=lambda i: i[1])]

                for name, angle in angs:
                    a = "{:3.1f}".format(angle)
                    self.add_row('', '', '', name, a)

                self.add_blank_row()
