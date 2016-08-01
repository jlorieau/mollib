"""
Creates reports for MolLib.
"""
from Report import ReportRenderer, SectionRenderer

class MoleculeReportRenderer(ReportRenderer):

    pass

class BasicStatisticsSection(SectionRenderer):

    def __init__(self, molecule):
        self.molecule = molecule
