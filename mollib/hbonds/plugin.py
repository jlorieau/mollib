"""
The plugin for the hbond submodule.
"""

from mollib.plugins import Plugin
from mollib.hbonds import find_hbond_partners, settings
from mollib.utils import MDTable

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

        if args.detailed:
            # Setup the table
            table = MDTable('Num', 'Donor', 'Acceptor',
                            'Parameter', 'Value')
            table.title = ('Hydrogen bond detailed '
                           'listing for {}'.format(molecule.name))

            # Add the Hbonds to the table
            for count, hbond in enumerate(hbonds, 1):
                # Get the dipole atom distances
                dists = [(atoms, d)
                         for atoms,d in sorted(hbond.distances.items(),
                                               key=lambda i: i[1])]

                # Convert the distance names and distances to strings
                # Add the atom names
                for i in range(len(dists)):
                    a = dists[i][0]
                    dist = dists[i][1]
                    name = ''.join(('{', a[0:2], '}...{', a[2:4], '}'))
                    name = name.format(a1=hbond.acceptor.atom1,
                                       a2=hbond.acceptor.atom2,
                                       d1=hbond.donor.atom1,
                                       d2=hbond.donor.atom2)

                    dist = '{:3.2f}'.format(dist)
                    dists[i] = (name, dist)

                # Get the dipole angles
                angs = [(name, a)
                         for name,a in sorted(hbond.angles.items(),
                                               key=lambda i: i[1])]

                # Add the hbond row with the first distance.
                table.add_row(count, hbond.donor, hbond.acceptor,
                              dists[0][0], dists[0][1])

                # Add the distance rows
                for name, d in dists[1:]:

                    # Create the rows
                    table.add_row('', '', '', name, d)

                # Add the angle rows
                for name, a in angs:
                    table.add_row('', '', '', name, a)

                table.add_blank_row()

            print(table.content())
        else:
            # Setup the table
            table = MDTable('Num', 'Donor', 'Acceptor', 'type (major)')
            table.title = ('Hydrogen bond '
                           'listing for {}'.format(molecule.name))

            # Add the Hbonds to the table
            for count, hbond in enumerate(hbonds, 1):
                table.add_row(count, hbond.donor, hbond.acceptor,
                              ' '.join((hbond.major_classification,
                                        hbond.minor_classification)))
            print(table.content())