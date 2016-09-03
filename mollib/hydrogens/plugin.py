"""
The plugin for the hydrogen submodule.
"""

import logging

from mollib.plugins import Plugin
from mollib.hydrogens import add_hydrogens


class Hydrogenate(Plugin):
    """The core plugin to offer the '--hydrogenate' preprocessor."""

    enabled = True
    order = -10  # process before the Process plugin

    def options(self, subparsers):
        # Get or create the parent 'process' parser
        parent = Plugin.parents.setdefault('process',
                                    subparsers.add_parser('', add_help=False))

        # Hydrogenate parameter
        parent.add_argument('--hydrogenate',
                             action='store_true', required=False,
                             help=("Strip hydrogens and re-add them before "
                                   "analysis"))

        return subparsers

    def preprocess(self, molecule, args):
        """Preprocess the molecule by adding hydrogens to it.
        """
        if 'hydrogenate' in args:
            add_hydrogens(molecule, strip=True)