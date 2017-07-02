# -*- coding: utf-8 -*-
"""
The plugin for the pa submodule.
"""
import os.path
from textwrap import TextWrapper

from mollib.plugins import Plugin
from mollib.utils.checks import check_not_empty
from mollib.utils.files import write_file
from mollib.utils.text import word_list
from mollib.utils.net import get_or_fetch
from mollib.utils.interactions import interaction_type
import mollib.utils.settings as utils_settings

from .data_readers import read_pa_file
from .process_molecule import Process
from .svd import calc_pa_SVD
from .fixers import Fixer
from .reports import report_tables, stats_table
from . import settings


# TODO: add '--exclude' option to remove couplings from data
# TDOO: add '--set' option to select RDC datasets
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

        p.add_argument('-d', '--data',
                       action='append', nargs='+',
                       required=True, metavar='id/filename',
                       help="(required) Alignment file or identifier with RDC "
                            "and RACS data")

        # Allow for the optional output of the results to a file
        p.add_argument('-o', '--out',
                       action='store', required=False,
                       type=str, metavar='filename',
                       help="The output filename for the reports of the fit "
                            "data.")

        p.add_argument('-p', '--pred',
                       action='store', required=False,
                       type=str, metavar='filename',
                       help='The output filename for the report of the ' 
                            'back-calculated RDCs and RACSs that are not in '
                            'the experimental data.')

        p.add_argument('--summary',
                       action='store_true',
                       required=False,
                       help="Only display the fit summary")

        p.add_argument('--set',
                       action='store', required=False,
                       metavar='id',
                       help='If multiple datasets are available, this option '
                            'specifies which dataset to use.')

        p.add_argument('--exclude',
                       action='store', required=False, nargs='*',
                       metavar='interaction-type',
                       help='Exclude one or more interactions of the following '
                            'type(s). ex: N-H or CE-HE')

        p.add_argument('--project-methyls',
                       action='store_true',
                       help='Fit methyl RDCs by projecting their values on the '
                            'corresponding C-C bond, as used by Xplor-NIH')

        p.add_argument('--methyl-scale',
                       action='store', required=False,
                       type=float, metavar='number',
                       help='The order parameter to use in scaling the methyl '
                            'RDCs.')

        fixers = p.add_argument_group("fixer arguments")
        fix_sign = fixers.add_mutually_exclusive_group()
        fix_sign.add_argument('--fix-sign',
                              action='store_true',
                              help="Check and fix mistakes in RDC and RACS "
                                   "sign")
        fix_sign.add_argument('--nofix-sign',
                              action='store_true',
                              help="Disable check in RDC and RACS sign")

        fix_nh_scale = fixers.add_mutually_exclusive_group()
        fix_nh_scale.add_argument('--fix-nh-scale',
                                  action='store_true',
                                  help="Check and rescale couplings that "
                                       "were scaled to the N-H RDC.")
        fix_nh_scale.add_argument('--nofix-nh-scale',
                                  action='store_true',
                                  help="Disable N-H rescaling of couplings.")

        # The following options can be turned off and on
        fix_outliers = fixers.add_mutually_exclusive_group()
        fix_outliers.add_argument('--fix-outliers',
                              action='store_true',
                              help="Fit without outliers")
        fix_outliers.add_argument('--nofix-outliers',
                              action='store_true',
                              help="Disable fitting without outliers")

    def help(self):
        return ("Singular Value Decomposition of partially aligned NMR data "
                "(RDCs and RACSs)")

    def preprocess(self, molecules, args):
        """Conduct argument checks."""
        # The --out and --pred filenames cannot be the same when using the 'pa'
        # subcommand
        if args.command == 'pa':
            if getattr(args, 'out', False) and args.out == args.pred:
                # Add '_pred' to the output filename
                args.pred = '_pred'.join(os.path.splitext(args.out))

    def process(self, molecules, args):
        """Process the SVD of molecules."""
        # Setup the configuration options
        if 'project_methyls' in args and args.project_methyls:
            settings.project_methyls = True
        if 'methyl_scale' in args and args.methyl_scale is not None:
            settings.methyl_order_parameter = args.methyl_scale
        if 'fix_sign' in args and args.fix_sign:
            settings.enable_signfixer = True
        if 'nofix_sign' in args and args.nofix_sign:
            settings.enable_signfixer = False
        if 'fix_nh_scale' in args and args.fix_nh_scale:
            settings.enable_nhscalefixer = True
        if 'nofix_nh_scale' in args and args.nofix_nh_scale:
            settings.enable_nhscalefixer = False
        if 'fix_outliers' in args and args.fix_outliers:
            settings.enable_outlierfixer = True
        if 'nofix_outliers' in args and args.nofix_outliers:
            settings.enable_outlierfixer = False

        # If specified, get the identifier for the dataset to use.
        set_id = args.set if 'set' in args else None

        # Process the partial alignment calculation
        if args.command == 'pa':
            # Get the alignment data
            data = {}
            for data_filename in args.data[0]:
                # verify that the file exists
                file_path = get_or_fetch(data_filename, extensions='mr.gz',
                                         urls=settings.mr_urls,
                                         critical=True)

                # Read the data from the file.
                data_dict = read_pa_file(file_path, set_id)
                data.update(data_dict)

            # If excluded interactions are specified, remove these.
            if args.exclude:
                data = {k:v for k, v in data.items()
                        if interaction_type(k) not in args.exclude}

            # verify that there is data in the data dict
            msg = "Could not find data in alignment data."
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

            # Prepare table of stats and fit values
            table = stats_table(stats)

            # Prepare a table of the observed and predicted data
            tables = report_tables(data, data_pred)

            if len(molecules) > 1:
                # Make title for stats table
                title = "Summary SVD Statistics for Molecules "
                title += word_list([m.fullname for m in molecules])
                table.title = title

                # Make title for the fit data table
                title = "Observed and Predicted RDCs and RACS for Molecules "
                title += word_list([m.fullname for m in molecules])
                tables['fit'].title = title

                # Make title for the back-calculated predicted data
                title = "Back-calculated RDCs and RACS for Molecules "
                title += word_list([m.fullname for m in molecules])
                tables['pred'].title = title
            else:
                # Make title for stats table
                title = "Summary SVD Statistics for Molecule "
                title += molecules[0].fullname
                table.title = title

                # Make title for the fit data table
                title = "Observed and Predicted RDCs and RACS for Molecule "
                title += molecules[0].fullname
                tables['fit'].title = title

                # Make title for the back-calculated predicted data
                title = "Back-calculated RDCs and RACS for Molecule "
                title += molecules[0].fullname
                tables['pred'].title = title

            # Prepare the standard output
            summary = table.content()
            output = tables['fit'].content()

            # Prepare and format the fixes listing
            if fixes:
                # Setup the text wrapper so that the lines of fixes do not
                # exceed the set maximum number of columns.
                wrapper = TextWrapper()
                wrapper.initial_indent = '* '
                wrapper.subsequent_indent = '  '
                wrapper.width = utils_settings.default_max_width

                fixes_wrapped = ['\n'.join(wrapper.wrap(fix)) for fix in fixes]
                fixes_output = '\n'.join(fixes_wrapped)
            else:
                fixes_output = ''

            # Print or write the report(s)
            print(summary)
            if args.out:
                output += fixes_output
                write_file('\n'.join((summary, output)), args.out)
            elif not args.summary:
                print(output)

            if fixes:
                print(fixes_output)

            # Write the predicted data
            if args.pred:
                write_file(tables['pred'].content(), args.pred)



