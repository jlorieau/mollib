"""
General Report Renderer

   @Author:             Justin L Lorieau <jlorieau>
   @Date:               2016-07-31T12:32:10-05:00
   @Last modified by:   jlorieau
   @Last modified time: 2016-08-01T09:28:27-05:00
   @License:            Copyright 2016
"""

import os
import subprocess
from section_renderers import *


class ReportRenderer(object):

    template_filename = 'template_header.md'
    output_filename = '{title}_report.pdf'

    title = None
    author = None

    def __init__(self, title, **kwargs):
        """Constructor

        [Required parameters]
            :title:     The title of the report.

        [Optional parameters]
            :template_filename: The filename of the template file used to render
                                the report.
            :output_filename:   The filename of the output file. If this is not
                                specified, it will be constructed from the
                                class's output_filename and title.
            :author:            Author of the report.
        """

        self.title = title

        if 'output_filename' not in kwargs:
            self.output_filename = self.output_filename \
                                       .format(title=title)
        else:
            self.output_filename = kwargs['output_filename']

        with open(self.template_filename) as f:
            self.template = f.read().format(title=title)

        for k,v in kwargs.values():
            if hasattr(self, k):
                setattr(self, k, v)

    def command_list(self):
        """Returns a list of commands and arguments to be executed by the
        renderer."""
        args = ['pandoc',]

        # Add variable definitions, if they're defined
        for attr in ['title', 'author']:
            value = getattr(self, attr, None)
            if value is not None and type(value) == str:
                args += ['-V', ':'.join((attr, value))]

        # Add the output argument
        args += ['-o', self.output_filename]
        return args

    def renderPDF(self, filename=None):
        filename = filename if filename is not None else self.output_filename

        args = self.command_list()

        process = subprocess.Popen(args=args, stdin=subprocess.PIPE)
        process.communicate(input=self.contents())

    def contents(self):
        contents = self.template

        return contents


### Tests ###
import unittest

class TestMolLib(unittest.TestCase):

    def test_report_renderer_header(self):
        report = ReportRenderer('2KXA')
        self.assertIn('2KXA', report.template)

    pass

    def test_render_pdf(self):
        report = ReportRenderer('2KXA')
        report.renderPDF()
        self.assertTrue(os.path.isfile(report.output_filename))

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
    unittest.main()
