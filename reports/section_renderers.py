"""
Section Renderer

   @Author:             Justin L Lorieau <jlorieau>
   @Date:               2016-08-01T08:01:00-05:00
   @Last modified by:   jlorieau
   @Last modified time: 2016-08-03T09:06:13-05:00
   @License:            Copyright 2016
"""
# imports for tests
import unittest


class SectionRenderer(object):
    """An Abstract Base Class for report sections.
    """

    def contents(self):
        """Returns the string contents of the section."""
        return ""


class HeaderSection(SectionRenderer):
    """Renders a section heading in Markdown.
    """

    def __init__(self, title, level=1):
        super(HeaderSection, self).__init__()
        self.title = title.strip()
        self.level = level

    def content(self):
        return ("\n" +
                "#" * self.level +
                " {title}".format(title=self.title) +
                "\n")


class TextSection(SectionRenderer):
    """Renders a basic markdown text section."""

    def __init__(self, text):
        super(TextSection, self).__init__()
        self.text = text

    def content(self):
        return self.text


class TableSection(SectionRenderer):
    """Renders a table section in Markdown.
    """

    def __init__(self, column_titles):
        """Constructor that initializes the column headers. The current
        implementation only uses single-line rows and headers. (No column
        width maximum)

        [Required Parameters]
            :column_titles:   An array of strings containing the column header
                              titles.
        """
        super(TableSection, self).__init__()
        self.column_titles = column_titles
        self.rows = []

    def add_row(self, row):
        """Add a row of values to the table"""
        # The number of items in the row must match the number of columns
        assert(len(row) == len(self.column_titles))

        # Add the rows and convert them to strings
        self.rows.append([str(i) for i in row])

    def content(self):
        """Renders a string for a markdown table.

        Note that this renders multiline tables in pandoc.
        """
        # Find the largest item for each column
        column_widths = [len(c) for c in self.column_titles]

        for col_index in range(len(column_widths)):
            max_width = column_widths[col_index]

            for row in self.rows:
                if len(row[col_index]) > max_width:
                    max_width = len(row[col_index])
            column_widths[col_index] = max_width

        # Format the table headers. All of the items are centered
        total_length = sum([width + 3 for width in column_widths]) - 1

        # Add top bar
        table = '\n' + '-' * total_length + '\n'

        # Add headers
        table += ' '.join([col.center(width + 2) for col, width in
                           zip(self.column_titles, column_widths)])
        table += '\n'

        # Add header bottom bars
        table += ''. join(['-' * (width + 2)
                           if count == 0 else ' ' + '-' * (width + 2)
                           for count, width in enumerate(column_widths)])
        table += '\n'

        # Add rows
        for row in self.rows:
            table += ' '.join([item.center(width + 2) for item, width in
                               zip(row, column_widths)])
            table += '\n\n'

        # Add bottom bar
        table += '-' * total_length + '\n\n'

        return table


class GraphSection(object):
    """Renders a graph section."""
    pass


# Tests


class TestSectionRenderer(unittest.TestCase):

    def test_header_section(self):
        "Tests the rendering of the HeaderSection."
        # Test level one header
        section = HeaderSection("My first section")
        self.assertEqual(section.content(),
                         "\n# My first section\n")

        # Test level two header. Test strip spaces.
        section = HeaderSection("     This header is level 2  ", level=2)
        self.assertEqual(section.content(),
                         "\n## This header is level 2\n")

    def test_text_section(self):
        "Test the rendering of the TextSection."
        text = """This is *my* text section with multiple
        lines and spaces."""

        section = TextSection(text)

        self.assertEqual(section.content(), text)

    def test_table_section(self):
        "Tests the rendering and spacing of the markdown table"
        section = TableSection(['one', 'two', 'three', 'four', 'and more'])
        section.add_row([1, 2, 3, 4, '+++'])
        section.add_row([1.0, 2.0, 3.0, 4.0, '+++'])
        section.add_row(['first item', 'second item', 'third item',
                         'a fourth item', 'an extra item'])
        content = section.content()

        target_content = ("\n"
"-----------------------------------------------------------------------\n"
"    one           two         three           four          and more   \n"
"------------ ------------- ------------ --------------- ---------------\n"
"     1             2            3              4              +++      \n\n"
"    1.0           2.0          3.0            4.0             +++      \n\n"
" first item   second item   third item   a fourth item   an extra item \n\n"
"-----------------------------------------------------------------------\n\n")

        self.assertEqual(content, target_content)