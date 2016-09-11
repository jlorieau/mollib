"""
The plugin for the hydrogen submodule.
"""

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

        # Hydrogenate parameter. Note that since the argument share the same
        # name as the class, this plugin will be active if the flag is passed.
        # (i.e. the selected method will return True)
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

    def selected(self, args):
        return True if args.hydrogenate else False