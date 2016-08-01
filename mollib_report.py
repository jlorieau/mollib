"""
Report Renderer for MolLib

   @Author:             Justin L Lorieau <jlorieau>
   @Date:               2016-07-31T12:32:10-05:00
   @Last modified by:   jlorieau
   @Last modified time: 2016-08-01T09:28:59-05:00
   @License:            Copyright 2016
"""


from Report import ReportRenderer, SectionRenderer

class MoleculeReportRenderer(ReportRenderer):

    pass

class BasicStatisticsSection(SectionRenderer):

    def __init__(self, molecule):
        self.molecule = molecule
