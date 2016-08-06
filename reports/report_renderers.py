"""
General Report Renderer

   @Author:             Justin L Lorieau <jlorieau>
   @Date:               2016-07-31T12:32:10-05:00
   @Last modified by:   jlorieau
   @Last modified time: 2016-08-02T21:36:38-05:00
   @License:            Copyright 2016
"""

import os
import subprocess
from .section_renderers import *

template = ("---\n"
            "title: {title}\n"
            "geometry: margin=1in\n"
            "---\n")

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
        self.renderers = []

        if 'output_filename' not in kwargs:
            self.output_filename = self.output_filename \
                                       .format(title=title)
        else:
            self.output_filename = kwargs['output_filename']

        self.template = template.format(title=title)

        for k,v in kwargs.values():
            if hasattr(self, k):
                setattr(self, k, v)

    def add_renderer(self, renderer):
        """Adds a Renderer to this Report Renderer."""
        self.renderers.append(renderer)

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
        process.communicate(input=self.content().encode())

    def content(self):
        content = self.template + '\n'
        content += '\n'.join([r.content() for r in self.renderers])
        content += '\n'
        return content


### Tests ###
import unittest

class TestMolLib(unittest.TestCase):

    def test_report_renderer_header(self):
        "Tests the proper placement of the title in the template."
        report = ReportRenderer('2KXA')
        self.assertIn('2KXA', report.template)

    def test_report_renderer_content(self):
        "Tests the integration of multiple SectionRenderers"
        report = ReportRenderer('My test report')

        # The report template should not depend on an external file to work.
        report.template = ('---\n'
                          'title: My test report\n'
                          'geometry: margin=1in\n'
                          '---\n')

        report.add_renderer(HeaderSection('My first title'))
        report.add_renderer(TextSection('This is *my* text\n\nHere!!'))
        report.add_renderer(HeaderSection('My sub-section', level=2))
        report.add_renderer(TextSection('With its own _text_.'))

        content = report.content()

        target_content = ('---\n'
                          'title: My test report\n'
                          'geometry: margin=1in\n'
                          '---\n'
                          '\n\n'
                          '# My first title\n\n'
                          'This is *my* text\n\n'
                          'Here!!\n\n'
                          '## My sub-section\n\n'
                          'With its own _text_.\n')

        self.assertEqual(content, target_content)


    def test_render_pdf(self):
        "Tests the create of the output pdf file"
        report = ReportRenderer('2KXA')
        report.renderPDF()
        self.assertTrue(os.path.isfile(report.output_filename))
