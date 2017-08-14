# -*- coding: utf-8 -*-
"""
The plugin for 'measure' command.
"""

from math import floor, log10
import itertools

import numpy as np

from mollib.plugins import Plugin
from mollib.core.geometry import (measure_distances, measure_angles,
                                  measure_dihedrals)
from mollib.utils import MDTable


def round_to_sigs(value, error):
    """Find the significant figures for error, and return a tuple with the
    value and error rounded to those sig figs."""
    # Determine the sig figs in the error so that the number can be properly
    # rounded
    try:
        sigs = -int(floor(log10(abs(error))))
    except ValueError:
        return (value, error)

    return (round(value, sigs), round(error, sigs))


# TODO: Implement the '--within' option.
class Measure(Plugin):
    """The core plugin to offer the 'measure' command."""

    enabled = True
    order = 100
    create_command_subparser = True

    def process_parser(self):
        "Process the parser for the 'measure' command."
        subparser = self.command_subparsers['measure']

        measure_group = subparser.add_argument_group('measurement options')

        # The following options can only be used one at a time (i.e. they're
        # mutually exclusive
        exc_group = measure_group.add_mutually_exclusive_group(required=False)

        exc_group.add_argument('-d', '--dist', nargs=2, required=False,
                               metavar='atom', type=str,
                               action='append',
                               help="Measure distances between 2 atom "
                                    "selections. ex: 31.N 32.CA")

        exc_group.add_argument('-a', '--angle', nargs=3, required=False,
                               metavar='atom', type=str,
                               action='append',
                               help="Measure angles between 3 atom selections. "
                                    "ex: 31.N 31.CA 31.C")

        exc_group.add_argument('-dih', '--dihedral', nargs=4, required=False,
                               metavar='atom', type=str,
                               action='append',
                               help="Measure dihedral angles between 4 atom "
                                    "selections. ex: 31.N 31.CA 31.C 32.N")

        exc_group.add_argument('-w', '--within', nargs=2, required=False,
                               metavar='atom', type=str,
                               action='append',
                               help="Measure all distances from atom selection "
                                    "to within the specified distance. "
                                    "ex: 31:33.N 5")

        measure_group.add_argument('--stats',
                                   required=False, action='store_true',
                                   help="Report statistics on the reported "
                                        "measurements.")

        # Arguments to filter the results
        filters = subparser.add_argument_group(title='filters')
        filters.add_argument('--only-intra',
                             dest='only_intra',
                             action='store_true', default=False,
                             help='Only report measurements within a residue')
        filters.add_argument('--exclude-intra',
                             dest='exclude_intra',
                             action='store_true', default=False,
                             help='Exclude measurements within a residue')
        filters.add_argument('--only-intra-chain',
                             dest='only_intra_chain',
                             action='store_true', default=False,
                             help='Only report measurements within a chain')
        filters.add_argument('--exclude-intra-chain',
                             dest='exclude_intra_chain',
                             action='store_true', default=False,
                             help='Exclude measurements within a chain')
        filters.add_argument('--only-delta',
                             dest='residue_delta', action='store',
                             default=None, metavar='DELTA', type=int,
                             help=('Only report residues separated by DELTA '
                                   'residue numbers'))
        filters.add_argument('--only-bonded',
                             dest='bonded',
                             action='store_true', default=None,
                             help=('Only report measurements from bonded '
                                   'atoms'))

    def help(self):
        return "Measure geometries in molecules"

    def process(self, molecules, args):
        """Measure geometries in molecules."""
        for molecule in molecules:
            self.process_molecule(molecule, args)

    def process_molecule(self, molecule, args):
        """Measure geometries in a given molecule."""
        if args.command == 'measure':
            if getattr(args, 'dist', False):
                table = MDTable('Num', 'Atom 1', 'Atom 2', 'Dist. (A)')
                table.title = 'Distances for {}'.format(molecule.fullname)

                dists_list = []
                for selector1, selector2 in args.dist:
                    # Get the distances
                    dists = measure_distances(molecule, selector1, selector2,
                                  only_intra=args.only_intra,
                                  exclude_intra=args.exclude_intra,
                                  only_intra_chain=args.only_intra_chain,
                                  exclude_intra_chain=args.exclude_intra_chain,
                                  residue_delta=args.residue_delta,
                                  bonded=args.bonded)
                    dists_list.append(dists)

                # Add the distances to the table
                for count, dist in enumerate(itertools.chain(*dists_list), 1):
                    atom1, atom2, d = dist
                    table.add_row(count, atom1, atom2, '{:.2f}'.format(d))

                # Calculate the stats if specified
                if args.stats:
                    values = [i[2] for i in itertools.chain(*dists_list)]
                    mean = np.mean(values)
                    std = np.std(values)
                    stat_str = '{} ± {}'.format(*round_to_sigs(mean, std))
                    stat_bar = '-' * len(stat_str)

                    table.add_row('', '', '', stat_bar)  # Draw a '----' bar
                    table.add_row('', '', '', stat_str)  # Add the stats

                # Print the table
                print(table.content())

            if getattr(args, 'angle', False):
                table = MDTable('Num', 'Atom 1', 'Atom 2', 'Atom 3',
                                'Angle (deg)')
                table.title = 'Angles for {}'.format(molecule.fullname)

                angs_list = []

                for selector1, selector2, selector3 in args.angle:
                    # get the angles
                    angs = measure_angles(molecule, selector1, selector2,
                                  selector3,
                                  only_intra=args.only_intra,
                                  exclude_intra=args.exclude_intra,
                                  only_intra_chain=args.only_intra_chain,
                                  exclude_intra_chain=args.exclude_intra_chain,
                                  residue_delta=args.residue_delta,
                                  bonded=args.bonded)
                    angs_list.append(angs)

                # Add the angles to the table
                for count, ang in enumerate(itertools.chain(*angs_list), 1):
                    atom1, atom2, atom3, a = ang
                    table.add_row(count, atom1, atom2, atom3,
                                  '{:.1f}'.format(a))

                # Calculate the stats if specified
                if args.stats:
                    values = [i[3] for i in itertools.chain(*angs_list)]
                    mean = np.mean(values)
                    std = np.std(values)
                    stat_str = '{} ± {}'.format(*round_to_sigs(mean, std))
                    stat_bar = '-' * len(stat_str)

                    table.add_row('', '', '', '', stat_bar)  # Draw a '----' bar
                    table.add_row('', '', '', '', stat_str)  # Add the stats

                # Print the table
                print(table.content())

            if getattr(args, 'dihedral', False):
                table = MDTable('Num', 'Atom 1', 'Atom 2', 'Atom 3', 'Atom 4',
                                'Dihedral (deg)')
                table.title = 'Dihedrals for {}'.format(molecule.fullname)

                dihs_list = []

                for selector1, selector2, selector3, selector4 in args.dihedral:
                    # Get the dihedrals
                    dihs = measure_dihedrals(molecule, selector1, selector2,
                                  selector3, selector4,
                                  only_intra=args.only_intra,
                                  exclude_intra=args.exclude_intra,
                                  only_intra_chain=args.only_intra_chain,
                                  exclude_intra_chain=args.exclude_intra_chain,
                                  residue_delta=args.residue_delta,
                                  bonded=args.bonded)
                    dihs_list.append(dihs)

                # Add the angles to the table
                for count, dih in enumerate(itertools.chain(*dihs_list), 1):
                    atom1, atom2, atom3, atom4, d = dih
                    table.add_row(count, atom1, atom2, atom3, atom4,
                                  '{:.1f}'.format(d))

                # Calculate the stats if specified
                if args.stats:
                    values = [i[4] for i in itertools.chain(*dihs_list)]
                    mean = np.mean(values)
                    std = np.std(values)
                    stat_str = '{} ± {}'.format(*round_to_sigs(mean, std))
                    stat_bar = '-' * len(stat_str)

                    table.add_row('', '', '', '', '',
                                  stat_bar)  # Draw a '----' bar
                    table.add_row('', '', '', '', '',
                                  stat_str)  # Add the stats

                # Print the table
                print(table.content())
