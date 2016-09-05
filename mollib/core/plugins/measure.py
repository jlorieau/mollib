"""
The plugin for 'measure' command.
"""

import logging

from mollib.plugins import Plugin
from mollib.core.geometry import (measure_distances, measure_angle,
                                  measure_dihedral)


class Measure(Plugin):
    """The core plugin to offer the 'measure' command."""

    enabled = True
    order = 100

    def options(self, subparsers):
        parser = super(Measure, self).options(subparsers)
        group = parser.add_mutually_exclusive_group(required=True)

        group.add_argument('-d', '--dist', nargs=2, required=False,
                        metavar='atom', type=str,
                        action='append',
                        help=("Measure distances between 2 atoms."
                              "ex: 31-N 32-CA"))

        group.add_argument('-a', '--angle', nargs=3, required=False,
                           metavar='atom', type=str,
                           action='append',
                           help=("Measure angles between 3 atoms."
                                 "ex: 31-N 31-CA 31-C"))

        group.add_argument('-w', '--within', nargs=2, required=False,
                           metavar='atom', type=str,
                           action='append',
                           help=("Measure all distances from atom to within "
                                 "the specified distance. "
                                 "ex: 31:33-N 5"))

        # Arguments to set inter-residue/intra-residue and chain options
        group2 = parser.add_mutually_exclusive_group(required=False)
        group2.add_argument('--exclude-intra',
                            dest='exclude_intra',
                            action='store_true', default=False,
                            help='Exclude intra-residue measurements')
        group2.add_argument('--delta',
                            dest='residue_delta', action='store',
                            default=None, metavar='DELTA', type=int,
                            help=('Only report residues separated by DELTA '
                                  'residue numbers.'))

        return parser

    def help(self):
        return "Measure geometries in molecules"

    def process(self, molecule, args):
        "Measure geometries in molecules."
        if 'dist' in args and args.dist is not None:
            msg = "({molecule}) {a1: >8} {a2: <8}: {dist:.2f} A"
            observed_pairs = {}

            for selector1, selector2 in args.dist:
                # Get the distances
                dists = measure_distances(molecule, selector1, selector2,
                                          residue_delta=args.residue_delta,
                                          exclude_intra=args.exclude_intra)

                # print the output message
                for dist in dists:
                    atom1, atom2, d = dist
                    print(msg.format(molecule=molecule.name,
                                     a1=atom1, a2=atom2, dist=d))

        if 'angle' in args and args.angle is not None:
            msg = "({molecule}) {a1: >8} {a2: ^8} {a3: <8}: {angle:.1f} deg"
            for a1, a2, a3 in args.angle:
                # This function logs an error if a1, a2 or a3 aren't properly
                # formatted.
                atoms1 = molecule.get_atoms(a1)
                atoms2 = molecule.get_atoms(a2)
                atoms3 = molecule.get_atoms(a3)

                if not atoms1 or not atoms2 or not atoms3:
                    continue

                for i in atoms1:
                   for j in atoms2:
                       for k  in atoms3:
                            # Skip if any atoms are the same
                            if i==j or j==k or i==k:
                                continue

                            # measure and print the angle
                            angle = measure_angle(i, j, k)
                            print(msg.format(molecule=molecule.name,
                                             a1=i, a2=j, a3=k,
                                             angle=angle))
