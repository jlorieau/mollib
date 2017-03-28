# -*- coding: utf-8 -*-
"""
The plugin for the hbond submodule.
"""
import os.path
import logging

from mollib.plugins import Plugin
from mollib.utils.checks import check_file, check_not_empty

from .data_readers import read_pa_file
from .process_molecule import Process
from .svd import calc_pa_SVD
from .reports import report_tables, stats_table


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
                check_file(data_filename, critical=True)

                # Read the data from the file
                data_dict = read_pa_file(data_filename)
                data.update(data_dict)

            # verify that there is data in the data dict
            msg = "Could not find data in alignment files."
            check_not_empty(data=data, msg=msg, critical=True)

            # Prepare the magnetic interactions for the molecules
            labels = data.keys()
            process = Process(molecules)
            magnetic_interactions = process.process(labels=labels)

            # Conduct the SVD on the data
            (data_pred, Saupe_components,
             stats) = calc_pa_SVD(magnetic_interactions, data)

            # Prepare a table of the observed and predicted data
            tables = report_tables(data, data_pred)
            # tables['fit'].title = "Molecule fit"

            # Prepare a table of the stats
            print(stats_table(stats).content())

            print(tables['fit'].content())

