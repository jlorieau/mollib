"""
The plugin for the hydrogen submodule.
"""

from mollib.plugins import Plugin
from mollib.hydrogens import add_hydrogens


class Hydrogenate(Plugin):
    """The core plugin to offer the '--hydrogenate' preprocessor."""

    enabled = True
    order = -10  # process before the Process plugin
    create_command_subparser = False

    def process_parser(self):
        "Add the '--hydrogenate' option to all subparsers."
        for subparser in self.command_subparsers.values():
            # Hydrogenate parameter. Note that since the argument shares the
            # same name as the class, this plugin will be active if the flag is
            # passed. (i.e. the selected method will return True)
            subparser.add_argument('--hydrogenate',
                                action='store_true', required=False,
                                help=("Strip hydrogens and re-add them before "
                                      "analysis"))

    def preprocess(self, molecules, args):
        """Preprocess molecules"""
        for molecule in molecules:
            self.preprocess_molecule(molecule, args)

    def preprocess_molecule(self, molecule, args):
        """Preprocess the molecule by adding hydrogens to it.
        """
        if getattr(args, 'hydrogenate', False):
            add_hydrogens(molecule, strip=True)

    def selected(self, args):
        return True if args.hydrogenate else False