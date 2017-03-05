# -*- coding: utf-8 -*-
"""
The plugin for the hbond submodule.
"""
import os.path
import logging

from mollib.plugins import Plugin

from .data_readers import read_pa_file
from .process_molecule import Process
from .svd import calc_pa_SVD
from .reports import report_tables


class PA(Plugin):
    """The core plugin to offer the 'pa' command."""

    enabled = True
    order = 110
    create_command_subparser = True

    def help(self):
        return "Fit partial alignment RDC and RACS data using SVD"

    def process_parser(self):
        "Process the parser for the 'pa' command."
        p = self.command_subparsers['pa']

        p.add_argument('-a', '--alignment',
                       action='append', nargs='+',
                       required=True,
                       help="Alignment file with RDC and RACS data")

        p.add_argument('-p', '--pred',
                       action='store_true',
                       help="Report predicted RDCs and RACS")

    def process(self, molecules, args):
        """Process the SVD of molecules."""
        if args.command == 'pa':
            # Get the alignment data
            data = {}
            for data_filename in args.alignment[0]:
                # verify that the file exists
                if not os.path.isfile(data_filename):
                    msg = "Filename '{}' does not exist."
                    logging.error(msg.format(data_filename))
                    continue

                # Read the data from the file
                data_dict = read_pa_file(data_filename)
                data.update(data_dict)

            # Prepare the magnetic interactions for the molecules
            process = Process(molecules)
            magnetic_interactions = process.process()

            # Conduct the SVD on the data
            (data_pred, Saupe_components,
             stats) = calc_pa_SVD(magnetic_interactions, data)

            # Report the statistics
            tables = report_tables(data, data_pred)
            # tables['fit'].title = "Molecule fit"

            print(tables['fit'].content())

