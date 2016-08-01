"""
General report renderer.

author: J Lorieau

Copyright 2016
"""
import os
import subprocess

class SectionRenderer(object):

    name = None
    header = "## {name}"

    def __init__(self, name):
        self.name = name
        self.header = self.header.format(name=name)

    def contents(self):
        return ""


class ReportRenderer(object):

    report_name = None
    header_filename = 'template_header.md'
    output_filename = '{report_name}_report.pdf'

    self.variables = {'author': None}


    def __init__(self, report_name, **kwargs):
        self.report_name = report_name
        self.output_filename = self.output_filename \
                                   .format(report_name=report_name)
        with open(self.header_filename) as f:
            self.header = f.read().format(report_name = report_name)

        for k,v in kwargs.values():
            if hasattr(self, k):
                setattr(self, k, v)

        self.variables = {}

        self.command = ["pandoc",]
        if self.author is not None:
            self.command += ["-V", "author:{author}".format(author=self.author)]

    def contents(self):
        contents = self.header

        return contents

    def renderPDF(self, filename=None):
        filename = filename if filename is not None else self.output_filename

        args = self.command + \
                ["-o", "{filename}".format(filename=filename)]
        process = subprocess.Popen(args=args, stdin=subprocess.PIPE)
        process.communicate(input=self.contents())

class PandocRenderer(ReportRenderer):

    @staticmethod
    def test_renderer():
        pass

### Tests ###
import unittest

class TestMolLib(unittest.TestCase):

    def test_report_renderer_header(self):
        report = ReportRenderer('2KXA')
        self.assertIn('2KXA', report.header)

    pass

    def test_render_pdf(self):
        report = ReportRenderer('2KXA')
        report.renderPDF()
        self.assertTrue(os.path.isfile(report.output_filename))

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
    unittest.main()
