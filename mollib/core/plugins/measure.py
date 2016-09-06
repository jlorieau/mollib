"""
The plugin for 'measure' command.
"""

import logging

from mollib.plugins import Plugin
from mollib.core.geometry import (measure_distances, measure_angles,
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

        # Arguments to filter the results
        #group2 = parser.add_mutually_exclusive_group(required=False)
        group2 = parser.add_argument_group(title='filters')
        group2.add_argument('--only-intra',
                            dest='only_intra',
                            action='store_true', default=False,
                            help='Only report measurements within a residue')
        group2.add_argument('--exclude-intra',
                            dest='exclude_intra',
                            action='store_true', default=False,
                            help='Exclude measurements within a residue')
        group2.add_argument('--only-intra-chain',
                            dest='only_intra_chain',
                            action='store_true', default=False,
                            help='Only report measurements within a chain')
        group2.add_argument('--exclude-intra-chain',
                            dest='exclude_intra_chain',
                            action='store_true', default=False,
                            help='Exclude measurements within a chain')
        group2.add_argument('--only-delta',
                            dest='residue_delta', action='store',
                            default=None, metavar='DELTA', type=int,
                            help=('Only report residues separated by DELTA '
                                  'residue numbers'))
        group2.add_argument('--only-bonded',
                            dest='bonded',
                            action='store_true', default=None,
                            help=('Only report measurements from bonded '
                                  'atoms'))

        return parser

    def help(self):
        return "Measure geometries in molecules"

    def process(self, molecule, args):
        "Measure geometries in molecules."
        # setup the filters

        if 'dist' in args and args.dist is not None:
            msg = "({molecule}) {a1: >8} {a2: <8}: {dist:.2f} A"

            for selector1, selector2 in args.dist:
                # Get the distances
                dists = measure_distances(molecule, selector1, selector2,
                                  only_intra=args.only_intra,
                                  exclude_intra=args.exclude_intra,
                                  only_intra_chain=args.only_intra_chain,
                                  exclude_intra_chain=args.exclude_intra_chain,
                                  residue_delta=args.residue_delta,
                                  bonded=args.bonded)

                # print the output message
                for dist in dists:
                    atom1, atom2, d = dist
                    print(msg.format(molecule=molecule.name,
                                     a1=atom1, a2=atom2, dist=d))

        if 'angle' in args and args.angle is not None:
            msg = "({molecule}) {a1: >8} {a2: ^8} {a3: <8}: {angle:.1f} deg"

            for selector1, selector2, selector3 in args.angle:
                # Get the residue number arguments. By default, it is intra
                # residue.
                exclude_intra = args.exclude_intra

                residue_delta = args.residue_delta
                residue_delta = (0 if (residue_delta is None and
                                       not exclude_intra)
                                 else residue_delta)

                angs = measure_angles(molecule, selector1, selector2, selector3,
                                  only_intra=args.only_intra,
                                  exclude_intra=args.exclude_intra,
                                  only_intra_chain=args.only_intra_chain,
                                  exclude_intra_chain=args.exclude_intra_chain,
                                  residue_delta=args.residue_delta,
                                  bonded=args.bonded)

                # print the output message
                for ang in angs:
                    atom1, atom2, atom3, a = ang
                    print(msg.format(molecule=molecule.name,
                                     a1=atom1, a2=atom2, a3=atom3,
                                     angle=a))


