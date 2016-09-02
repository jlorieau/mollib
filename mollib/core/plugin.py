"""
The plugin for the core submodule.
"""

import argparse

from mollib.plugins import Plugin, check_number_arguments


class Process(Plugin):
    """The core plugin to offer the 'process' command."""

    enabled = True
    order = 0

    def options(self, subparsers):
        # Create the parent processor parser to be used by other plugins.
        parent = subparsers.add_parser('', add_help=False)
        parent.order = self.order

        # Input filename or identifier
        parent.add_argument('-i',
                            action='append', nargs='+', required=True, type=str,
                            metavar='id/filename',
                            help=("(required) The filename(s) or PDB "
                             "identifier(s) of the structure(s)"))

        # Output filename
        parent.add_argument('-o', '--out',
                            nargs=1, required=False, type=str,
                            metavar='filename',
                            help="The output file's name for the structure")

        # Config filename
        parent.add_argument('-c', '--config',
                            nargs=1, required=False, type=str,
                            metavar='filename',
                            help="The configuration file's name.")

        # Add it to the list of parent parsers. This has to be added to the
        # parent Plugin ABC to be seen by Plugin subclasses Sorting is needed
        # so that the order of options in the command line help message do not
        # change order.
        Plugin.parents = sorted(Plugin.parents + [parent,],
                                key=lambda p: p.order)

        # Create the process parser
        parser = subparsers.add_parser(self.command, parents=self.parents,
                                       help="Process a molecular structure")
        parser._optionals.title = self.argument_title

        return parser

    def help(self):
        return "Process the structure"

    def process(self, molecule):
        "Process and write the molecule file."
        pass

    def selected(self, args):
        "This plugin is always active."
        return True

class Measure(Plugin):
    """The core plugin to offer the 'measure' command."""

    enabled = True
    order = 10

    def options(self, subparsers):
        parser = super(Measure, self).options(subparsers)
        group = parser.add_mutually_exclusive_group(required=True)

        # One of the following has to be used
        group.add_argument('-a', '--atoms', nargs='+', required=False,
                        metavar='atom', type=str,
                        action=check_number_arguments(2, 4, 'atoms'),
                        help=("Measure distance, angle or dihedral between "
                              "atoms. 2 atoms gives distance, 3 atoms gives "
                              "angle, 4 atoms gives dihedral. "
                              "ex: N31-N C32-CA"))
        group.add_argument('-r', '--ramachandran', action='store_true',
                           required=False,
                           help=("Report a table of the Ramachandran angles."))
        return parser

    def help(self):
        return "Measure geometries in molecules"

    def process(self, molecule):
        "Measure geometries in molecules."
        pass