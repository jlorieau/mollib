"""
The plugin for the 'process' command.
"""

import logging

from mollib.plugins import Plugin


class Process(Plugin):
    """The core plugin to offer the 'process' command."""

    enabled = True
    order = 0
    create_command_subparser = True

    def process_parser(self):
        "Process the 'process' command parser."
        # Add the following commands to all subparsers
        for subparser in self.command_subparsers.values():
            # Input filename or identifier
            subparser.add_argument('-i', '--in', dest='i',
                                action='append', nargs='+', required=True,
                                type=str,
                                metavar='id/filename',
                                help=("(required) The filename(s) or PDB "
                                      "identifier(s) of the structure(s)"))

            # Output filename
            subparser.add_argument('-o', '--out',
                                action='append', nargs='*', required=False,
                                type=str, metavar='filename',
                                help="The output filename(s) for the structure(s)")

            # Config filename
            subparser.add_argument('-c', '--config',
                                nargs=1, required=False, type=str,
                                metavar='filename',
                                help="The configuration filename")

            # List molecule details
            subparser.add_argument('-l', '--list',
                                action='store_true',
                                help='List details on the molecule(s)')

    def help(self):
        return "Process the structure"

    def process(self, molecule, args):
        """Process the molecule

        - list the molecule details


        .. note:: Opening of the molecule (-i) and the configuration (-c) is
                  conducted elsewhere. Also, the file is written as the last
                  operation.
        """
        if args.list:
            print(molecule)
            chain_msg = '\tChain {:<3}: {:>4} residues, {:>4} atoms.'
            for chain in molecule.chains:
                print(chain_msg.format(chain, chain.residue_size,
                                       chain.atom_size))

    def postprocess(self, molecule, args):
        """Postprocess the molecules.

        - Writes the PDB file, if specified
        """
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
