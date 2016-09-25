# -*- coding: utf-8 -*-
"""
The plugin for the hbond submodule.
"""
from math import exp

import numpy as np

import mollib.core.settings
from mollib.plugins import Plugin
from mollib.hbonds import find_hbond_partners, classify_residues, settings
from mollib.utils import MDTable

class Hbonds(Plugin):
    """The core plugin to offer the 'Hbonds' command."""

    enabled = True
    order = 100
    create_command_subparser = True

    def help(self):
        return "Find and report hydrogen bonds in molecules"

    def process_parser(self):
        "Process the parser for the 'hbonds' command."
        p = self.command_subparsers['hbonds']

        p.add_argument('--aliphatic',
                       action='store_true',
                       help="Includes aliphatic hydrogen bonds")
        p.add_argument('--detailed',
                       action='store_true',
                       help="Report detailed information on hydrogen bonds.")
        p.add_argument('--sort-type',
                       action='store_true',
                       help='Sort hydrogen bonds by type')

        p = self.command_subparsers['measure']

        p.add_argument('--rama',
                       action='store_true',
                       help=("Report the Ramachandran angles. Filters and "
                             "options are ignored."))

    def process(self, molecule, args):
        """Process the molecule by finding and reporting its hydrogen bonds.
        """
        if args.command == 'hbonds':
            if getattr(args, 'aliphatic', False):
                settings.donor2_elements += "|C|13C"
                settings.hbond_distance_cutoff['d1a1']= (1.8, 3.0)

            hbonds = find_hbond_partners(molecule)

            # Sort the hbonds by the given criteria
            if getattr(args, 'sort_type', False):
                hbonds = sorted(hbonds,
                                key=lambda hb: (hb.major_classification,
                                                hb.minor_classification,
                                                hb.donor.atom2.residue.number))

            if getattr(args, 'detailed', False):
                # Setup the table
                table = MDTable('Num', 'Donor', 'Acceptor',
                                'Parameter', 'Value')
                table.title = ('Hydrogen bond detailed '
                               'listing for {}'.format(molecule.name))

                # Add the Hbonds to the table
                for count, hbond in enumerate(hbonds, 1):
                    # Get the dipole atom distances
                    dists = [(atoms, d)
                             for atoms,d in sorted(hbond.distances.items(),
                                                   key=lambda i: i[1])]

                    # Convert the distance names and distances to strings
                    # Add the atom names
                    for i in range(len(dists)):
                        a = dists[i][0]
                        dist = dists[i][1]
                        name = ''.join(('{', a[0:2], '}...{', a[2:4], '}'))
                        name = name.format(a1=hbond.acceptor.atom1,
                                           a2=hbond.acceptor.atom2,
                                           d1=hbond.donor.atom1,
                                           d2=hbond.donor.atom2)

                        dist = '{:3.2f}'.format(dist)
                        dists[i] = (name, dist)

                    # Get the dipole angles
                    angs = [(name, a)
                             for name,a in sorted(hbond.angles.items(),
                                                   key=lambda i: i[1])]

                    # Add the hbond row with the first distance.
                    table.add_row(count, hbond.donor, hbond.acceptor,
                                  dists[0][0], dists[0][1])

                    # Add the distance rows
                    for name, d in dists[1:]:

                        # Create the rows
                        table.add_row('', '', '', name, d)

                    # Add the angle rows
                    for name, a in angs:
                        table.add_row('', '', '', name, a)

                    table.add_blank_row()

                print(table.content())
            else:
                # Setup the table
                table = MDTable('Num', 'Donor', 'Acceptor', 'Classification')
                table.title = ('Hydrogen bond '
                               'listing for {}'.format(molecule.name))

                # Add the Hbonds to the table
                for count, hbond in enumerate(hbonds, 1):
                    table.add_row(count, hbond.donor, hbond.acceptor,
                                  '{} ({})'.format(hbond.major_classification,
                                                   hbond.minor_classification))
                print(table.content())

        # Process the Ramachandran angles. This function detects secondary
        # structure units from hbonds.
        if getattr(args, 'rama', False):
            # Setup the table
            table = MDTable('Residue', 'Phi (deg)', 'Psi (deg)',
                            'Classification', 'E (kT)/Prob.')
            table.title = ('Ramachandran angles '
                           'for {}'.format(molecule.name))

            # Classify the residues based on their backbone amide hydrogen
            # bonds
            classify_residues(molecule)

            # Populate the table with the ramachandran angles and secondary
            # structure classifications.
            energies = []
            for residue in molecule.residues:
                # Skip heteroatom chains
                if '*' in residue.chain.id:
                    continue

                phi, psi = residue.ramachandran_angles
                classification = getattr(residue, 'hbond_classification', '')
                energy = getattr(residue, 'energy_ramachandran', '-')

                # If the energy has a value (float) and it is above the energy
                # cutoff, add its value to the table. Otherwise, just print
                # a '-' character, if it is within acceptable ranges.
                if (isinstance(energy, float) and
                    energy > mollib.core.settings.energy_cutoff_ramachandran):
                    energies.append(energy)

                    # Convert to string.
                    # According to Morris. Proteins 12:345 (1992).
                    # Core: 81.9% population. E(kT) = 1.71
                    # Allowed: 14.8%. E(kT) = -3.41
                    # Generous: 2.0%. E(kT) = -4.34
                    E_prob = "{:>2.1f}/{:<3.1f}%"
                    E_prob = E_prob.format(energy, exp(-1. * energy) * 100.)
                else:
                    E_prob = '✓'.center(10)

                table.add_row('{}.{}'.format(residue.chain.id, residue),
                              "{:>6.1f}".format(phi or 0.),
                              "{:>6.1f}".format(psi or 0.),
                              classification,
                              E_prob)
            # energy_mean = np.mean(energies)
            # energy_std = np.std(energies)
            # table.add_row('', '', '', '',
            #               "{:>2.1f} +/- {:<2.1f}".format(energy_mean,
            #                                              energy_std))


            print(table.content())
