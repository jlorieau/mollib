"""
The plugin for the core submodule.
"""

import logging

from mollib.plugins import Plugin
from .geometry import measure_distance, measure_angle, measure_dihedral


class Process(Plugin):
    """The core plugin to offer the 'process' command."""

    enabled = True
    order = 0

    def options(self, subparsers):

        # Create the parent processor parser to be used by other plugins.
        parent = Plugin.parents.setdefault('process',
                                    subparsers.add_parser('', add_help=False))

        # Input filename or identifier
        parent.add_argument('-i', '--in', dest='i',
                            action='append', nargs='+', required=True, type=str,
                            metavar='id/filename',
                            help=("(required) The filename(s) or PDB "
                             "identifier(s) of the structure(s)"))

        # Output filename
        parent.add_argument('-o', '--out',
                            action='append', nargs='*', required=False,
                            type=str, metavar='filename',
                            help="The output filename(s) for the structure(s)")

        # Config filename
        parent.add_argument('-c', '--config',
                            nargs=1, required=False, type=str,
                            metavar='filename',
                            help="The configuration filename")


        # Create the process parser
        parser = subparsers.add_parser(self.command,
                                       parents=[parent,],
                                       help="Process a molecular structure")
        parser._optionals.title = self.argument_title

        return parser

    def help(self):
        return "Process the structure"

    def process(self, molecule, args):
        """Do nothing.

        Opening of the molecule (-i) and the configuration (-c) is
        conducted elsewhere. Also, the file is written as the last operation.
        """
        pass

    def postprocess(self, molecule, args):
        """Writes the PDB file, if specified"""
        # Do nothing if no output filename was given
        if args.out is None:
            return None

        # Get the corresponding filename
        output_filename = [o for i,o in zip(args.i[0], args.out[0])
                           if i==molecule.identifier]

        if len(output_filename) < 1:
            msg = "No output filename was specificed for {}."
            logging.error(msg.format(molecule.identifier))
            return None

        # Write the file.
        output_filename = output_filename[0]
        msg = "Writing ({}) to {}."
        logging.debug(msg.format(molecule.name, output_filename))

        molecule.write_pdb(output_filename)

    def selected(self, args):
        "This plugin is always active."
        return True

# def check_interresidue(delta=0, *args):
#     """Return try if the atoms listed in args are all within 'delta' number
#     of residues from each other.
#
#
#     """



class Measure(Plugin):
    """The core plugin to offer the 'measure' command."""

    enabled = True
    order = 100

    def options(self, subparsers):
        parser = super(Measure, self).options(subparsers)
        group = parser.add_mutually_exclusive_group(required=True)

        ## mutually exclusive arguments
        group.add_argument('-l', '--list',
                           action='store_true',
                           help="List details on the molecule")

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

        # group.add_argument('-r', '--ramachandran', action='store_true',
        #                    required=False,
        #                    help=("Report a table of the Ramachandran angles."))

        # Arguments avaliable to all other arguments
        parser.add_argument('--intra',
                           action='store_true',
                           help='Only report intraresidue measurements')

        parser.add_argument('--inter',
                            action='store_true',
                            help='Only report interresidue measurements')
        return parser

    def help(self):
        return "Measure geometries in molecules"

    def process(self, molecule, args):
        "Measure geometries in molecules."
        if 'dist' in args and args.dist is not None:
            msg = "({molecule}) {a1: >8} {a2: <8}: {dist:.2f} A"
            observed_pairs = {}

            for a1, a2 in args.dist:
                # This function logs an error if a1 or a2 isn't properly
                # formatted. An additional message is not needed. Just skip
                # it if both atoms aren't found.
                atoms1 = molecule.get_atoms(a1)
                atoms2 = molecule.get_atoms(a2)

                if not atoms1 or not atoms2:
                    continue

                for i in atoms1:
                    for j in atoms2:
                        # Skip if they're the same atom or if the distance
                        # has already be printed
                        if (i==j or
                            (observed_pairs.get(i, None) and
                             j in observed_pairs[i])):
                            continue

                        # Skip inter-residue if specified
                        if (args.intra and
                            (i.chain.id != j.chain.id or
                             i.residue.number != j.residue.number)):
                            continue

                        # Skip intra-residue if specified
                        if (args.inter and
                            (i.chain.id == j.chain.id and
                             i.residue.number == j.residue.number)):
                            continue

                        # Mark the pair as observed (to avoid duplicates)
                        observed_pairs.setdefault(j, []).append(i)

                        # measure and print the output message
                        dist = measure_distance(i, j)
                        print(msg.format(molecule=molecule.name, a1=i, a2=j,
                                         dist=dist))

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
