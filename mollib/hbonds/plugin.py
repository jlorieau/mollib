"""
The plugin for the hbond submodule.
"""

from mollib.plugins import Plugin
from mollib.hbonds import find_hbond_partners, settings


class Hbonds(Plugin):
    """The core plugin to offer the 'Hbonds' command."""

    enabled = True
    order = 100

    def help(self):
        return "Find and report hydrogen bonds in molecules"

    def options(self, subparsers):
        "Setup the argument parser"
        p = super(Hbonds, self).options(subparsers)

        p.add_argument('--aliphatic',
                       action='store_true',
                       help="Includes aliphatic hydrogen bonds")
        p.add_argument('--detailed',
                       action='store_true',
                       help="Report detailed information on hydrogen bonds.")
        return p

    def process(self, molecule, args):
        """Process the molecule by finding and reporting its hydrogen bonds.
        """
        if args.aliphatic:
            settings.donor2_elements += "|C|13C"
            settings.hbond_distance_cutoff['d1a1']= (1.8, 3.0)

        hbonds = find_hbond_partners(molecule)
        for count, hbond in enumerate(hbonds, 1):
            if args.detailed:
                print("{}. ({}) {}".format(count, molecule.name, hbond))
            else:
                print("{}. ({}) {}".format(count, molecule.name,
                                           hbond.short_repr()))
