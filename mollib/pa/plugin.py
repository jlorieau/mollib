# -*- coding: utf-8 -*-
"""
The plugin for the hbond submodule.
"""
import os.path
import logging

from mollib.plugins import Plugin
from mollib.utils.checks import check_file, check_not_empty
from mollib.utils.files import write_file
from mollib.utils.text import word_list

from .data_readers import read_pa_file
from .process_molecule import Process
from .svd import calc_pa_SVD
from .fixers import Fixer
from .reports import report_tables, stats_table
from . import settings


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
                       help="(required) Alignment file with RDC and RACS data")

        # Allow for the optional output of the results to a file
        p.add_argument('-o', '--out',
                       action='store', required=False,
                       type=str, metavar='filename',
                       help="The output filename for the reports")

        # The following options can be turned off and on
        fix_sign = p.add_mutually_exclusive_group()
        fix_sign.add_argument('--fix-sign',
                              action='store_true',
                              help="Check and fix mistakes in RDC and RACS "
                                   "sign")
        fix_sign.add_argument('--nofix-sign',
                              action='store_true',
                              help="Disable check in RDC and RACS sign")

        # The following options can be turned off and on
        fix_outliers = p.add_mutually_exclusive_group()
        fix_outliers.add_argument('--fix-outliers',
                              action='store_true',
                              help="Fit without outliers")
        fix_outliers.add_argument('--nofix-outliers',
                              action='store_true',
                              help="Disable fitting without outliers")


    def process(self, molecules, args):
        """Process the SVD of molecules."""
        # Setup the configuration options
        if 'fix_sign' in args and args.fix_sign:
            settings.enable_signfixer = True
        if 'nofix_sign' in args and args.nofix_sign:
            settings.enable_signfixer = False
        if 'fix_outliers' in args and args.fix_outliers:
            settings.enable_outlierfixer = True
        if 'nofix_outliers' in args and args.nofix_outliers:
            settings.enable_outlierfixer = False

        # Process the partial alignment calculation
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

            # Apply the fixers to see if the input data can be improved
            fixer = Fixer(molecules)
            data_fixed, fixes = fixer.fix(data)
            data = data_fixed if data_fixed is not None else data

            # Conduct the SVD on the data
            (data_pred, Saupe_components,
             stats) = calc_pa_SVD(magnetic_interactions, data)

            # Print the table of stats and fit values
            table = stats_table(stats)

            # Prepare a table of the observed and predicted data
            tables = report_tables(data, data_pred)

            if len(molecules) > 1:
                # Make title for stats table
                title = "Summary SVD Statistics for Molecules "
                title += word_list([m.name for m in molecules])
                table.title = title

                # Make title for pred and calc table
                title = "Observed and Predicted RDCs and RACS for Molecules "
                title += word_list([m.name for m in molecules])
                tables['fit'].title = title
            else:
                # Make title for stats table
                title = "Summary SVD Statistics for Molecule "
                title += molecules[0].name
                table.title = title

                # Make title for pred and calc table
                title = "Observed and Predicted RDCs and RACS for Molecule "
                title += molecules[0].name
                tables['fit'].title = title

            # Prepare the standard output
            output = '\n'.join((table.content(),
                                tables['fit'].content(),)
                               )
            if fixes:
                output += '\n'
                output += '\n'.join(['* ' + fix for fix in fixes])

            # Print or write the report(s)
            if args.out:
                write_file(output, args.out)
            else:
                print(output)

