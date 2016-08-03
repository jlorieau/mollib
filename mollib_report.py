"""
Report Renderer for MolLib

   @Author:             Justin L Lorieau <jlorieau>
   @Date:               2016-07-31T12:32:10-05:00
   @Last modified by:   jlorieau
   @Last modified time: 2016-08-02T21:37:10-05:00
   @License:            Copyright 2016
"""

import os
from .mollib import Molecule
from .Report import ReportRenderer
from .Report import HeaderSection, TextSection, TableSection


class MoleculeReportRenderer(ReportRenderer):
    "A Report for a MolLib molecule"

    def __init__(self, molecule, **kwargs):
        """Contructor for the Molecule report.

        [Required Parameters]
            :molecule:    A MolLib molecule to create the report for.
        """
        # The molecule name may contain path information like:
        # ../../../2KXA.pdb. This operation will strip the path information.
        name = os.path.split(molecule.name)[-1]

        # Construct a report title from the molecule's name.
        title = 'Report for molecule ' + name
        super(MoleculeReportRenderer,self).__init__(title=title, **kwargs)

        self.molecule = molecule
        self.build_report()

    def build_report(self):
        """Builds up the sections of the report."""

        self.add_renderer(HeaderSection("Protein Summary"))

        # Make a table on the masses of each chain
        columns = ['chain', '# of residues', 'mass (Da)']
        table = TableSection(column_titles=columns)
        for chain in self.molecule.chains:
            table.add_row([chain.id, chain.residue_size, chain.mass])
        self.add_renderer(table)


### Tests ###

import unittest

class TestMoleculeReportRenderer(unittest.TestCase):

    def test_basic_report(self):
        mol = Molecule('3C9J')
        report = MoleculeReportRenderer(molecule=mol)

        report.renderPDF()
