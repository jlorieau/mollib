# -*- coding: utf-8 -*-
"""
The plugin for the hbond submodule.
"""
from math import exp

import mollib.core.settings
from mollib.plugins import Plugin
from mollib.hbonds import find_hbond_partners, classify_residues, settings
from mollib.utils import MDTable, FormattedStr
from hbond_table import HBondTable


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
            # Change the hydrogen bond search patterns for aliphatic
            # hydrogen bonds
            if getattr(args, 'aliphatic', False):
                settings.donor2_elements += "|C|13C"
                settings.hbond_distance_cutoff['d1a1']= (1.8, 3.0)

            # Get the specified settings
            if hasattr(args, 'sort_type'):
                settings.hbond_table_sort_type = args.sort_type
            if hasattr(args, 'detailed'):
                settings.hbond_table_detailed = args.detailed

            # Measure the hydrogen bonds
            hbonds = find_hbond_partners(molecule)

            # Setup and print the table
            table = HBondTable(hbonds)
            table.title = ('Hydrogen bond listing for {}'.format(molecule.name))
            print(table.content())

        # Process the Ramachandran angles. This function detects secondary
        # structure units from hbonds.
        if getattr(args, 'rama', False):
            # Setup the table
            table = MDTable('Residue', 'Phi (deg)', 'Psi (deg)',
                            'Classification', 'E (kT) / Prob.')
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
                if isinstance(energy, float):
                    if energy < mollib.core.settings.energy_cutoff_good:
                        prob = exp(-1. * energy) * 100.
                        E_prob = "{:>2.1f} / {:>4.1f}%".format(energy, prob)
                        E_prob = FormattedStr(E_prob, 'green')
                    elif energy < mollib.core.settings.energy_cutoff_warning:
                        prob = exp(-1. * energy) * 100.
                        E_prob = "{:>2.1f} / {:>4.1f}%".format(energy, prob)
                        E_prob = FormattedStr(E_prob, 'yellow')

                    else:
                        prob = exp(-1. * energy) * 100.
                        E_prob = "{:>2.1f} / {:>4.1f}%".format(energy, prob)
                        E_prob = FormattedStr(E_prob, 'red')

                table.add_row('{}.{}'.format(residue.chain.id, residue),
                              "{:>6.1f}".format(phi or 0.),
                              "{:>6.1f}".format(psi or 0.),
                              classification,
                              E_prob)

            print(table.content())
