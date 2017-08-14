# -*- coding: utf-8 -*-
"""
The plugin for the hbond submodule.
"""
from mollib.plugins import Plugin
from mollib.hbonds import find_hbond_partners, settings
from .hbond_table import HBondTable
from .rama_table import RamaTable


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

        group = p.add_argument_group('hbond options')

        group.add_argument('--aliphatic',
                          action='store_true',
                           help="Includes aliphatic hydrogen bonds")
        group.add_argument('--detailed',
                           action='store_true',
                           help="Report detailed information on hydrogen bonds.")
        group.add_argument('--sort-type',
                           action='store_true',
                           help='Sort hydrogen bonds by type')

        p = self.command_subparsers['measure']

        # Find the 'measure' group and add the --rama option
        group = None
        for ag in p._action_groups:
            if ag.title.startswith('measure'):
                group = ag
        p = group if group is not None else p

        p.add_argument('--rama',
                       action='store_true',
                       help=("Report the Ramachandran angles. Filters and "
                             "options are ignored."))

    def process(self, molecules, args):
        """Process the molecule."""
        for molecule in molecules:
            self.process_molecule(molecule, args)

    def process_molecule(self, molecule, args):
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
            table.title = ('Hydrogen bond listing '
                           'for {}'.format(molecule.fullname))
            print(table.content())

        # Process the Ramachandran angles. This function detects secondary
        # structure units from hbonds.
        if getattr(args, 'rama', False):
            # Setup the table
            table = RamaTable(molecule)
            table.title = ('Ramachandran angles '
                           'for {}'.format(molecule.fullname))
            print(table.content())
