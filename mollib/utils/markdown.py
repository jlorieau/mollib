"""
Utilities for rendering information in Markdown
"""

import textwrap
from math import ceil

from . import settings


def print_lines(text_items, widths):
    """Given the unwrapped text_items, print them to fit in the given widths
    using textwrap.

    Parameters
    ----------
    text_items: list of str
        Items to print over multiple lines.
    widths: list of int
        width to fit the lines to.

    Examples
    --------
    >>> t = print_lines(['My', 'first and ever', 'multiline text'], \
                        [15, 15, 15])
    >>> print(t)
           My       first and ever multiline text
    >>> t = print_lines(['My', 'first and ever', 'multiline text'], \
                        [10, 10, 10])
    >>> t.splitlines()
    ['    My    first and multiline ', '             ever      text   ']
    """
    assert(len(text_items) == len(widths))

    # wrap the text_items
    text_items = [textwrap.wrap(t, w) for t,w in zip(text_items, widths)]
    no_lines = [len(t) for t in text_items]
    max_no_lines = max(no_lines)
    text = ''
    for line_no in range(max_no_lines):
        for text_item, width in zip(text_items, widths):
            text += (text_item[line_no].center(width)
                     if line_no < len(text_item) else
                     ''.center(width))

        # Print a new line on all lines except the last one
        if line_no != max_no_lines - 1 :
            text += '\n'
    return text


class MDTable(object):
    """Renders a table in Markdown.

    Parameters
    ----------
    title: str, optional
        If specified, this attribute will be used as the table's title.
    column_titles: iterable
        An iterable (list) of strings for the column header titles
    rows: iterable
        An iterable (list) of rows. Each row is an iterable (list) of strings
        that match the number of columns in column_titles
    multiline: bool
        If True, this table will render as a multiline table
    max_width: int, optional
        The maximum width of the table. If the table columns exceed this value,
        the table is converted to a multiline and columns are wrapped to
        accomodate the max_width.
    """

    title = None
    column_titles = None
    rows = None
    multiline = False
    max_width = settings.default_table_max_width

    def __init__(self, *args):
        """Constructor that initializes the column headers.

        The current implementation only uses single-line rows and headers.
        (No column columnwidth maximum)

        Parameters
        ----------
        *args: iterable
           A iterable (list) of strings for the column header titles
        """
        self.column_titles = args
        self.rows = []
        self.title = None


    def num_cols(self):
        "Return the number of columns in the table."
        return len(self.column_titles)

    def add_row(self, *args):
        """Add a row of values to the table"""
        # The number of items in the row must match the number of columns
        assert(len(args) == len(self.column_titles))

        # Add the rows and convert them to strings
        self.rows.append([str(i) for i in args])

    def add_blank_row(self):
        """Add a blank row to the table"""
        num_cols = self.num_cols()
        self.rows.append(['' for i in range(num_cols)])

    def column_widths(self):
        """Return the list of column widths for the current data.

        Returns
        -------
        list of tuples
            Each item in the returned list is a (width, wrap) tuple of the
            column_width (width, int) and whether to textwrap the column (wrap,
            bool)
        """
        # Find the largest item for each column
        column_widths = [len(c) for c in self.column_titles]

        for col_index in range(len(column_widths)):
            max_width = column_widths[col_index]

            for row in self.rows:
                if len(row[col_index]) > max_width:
                    max_width = len(row[col_index])
            column_widths[col_index] = max_width + 2  # Add a space around each
                                                      # Column

        # Now if the sum of the column_widths is larger than the self.max_width,
        # wrap the larger columns to fit the width
        total_width = sum(column_widths)
        if total_width > self.max_width:
            self.multiline = True
            difference = total_width - self.max_width
            difference_per_column = int(ceil(float(difference) /
                                             float(len(column_widths)))) + 2

            # Find which columns need to be textwrapped. The
            # largest_column_widths contains tuples of the index number(i) of
            # the corresponding column in column_widths and the new widths.
            largest_column_widths = [(i, w) for i,w in enumerate(column_widths)]
            largest_column_widths = sorted(largest_column_widths,
                                           key=lambda i: i[1], reverse=True)

            while(sum([c[1] for c in largest_column_widths]) > self.max_width):
                # Shrink the largest column by difference_per_column
                index, width = largest_column_widths[0]
                new_width = width - difference_per_column

                # Make sure the column isn't too small.
                if new_width < 5:
                    new_width = 5

                largest_column_widths[0] = (index, new_width)

                # Resort the largest_column_widths by index number, and copy over
                # to the column widths
                largest_column_widths = sorted(largest_column_widths,
                                               key= lambda i: i[1],
                                               reverse=True)

            for index, width in largest_column_widths:
                if width < column_widths[index]:
                    column_widths[index] = width

        return column_widths

    def content(self):
        """Renders a string for a markdown table.

        Note that this renders multiline tables in pandoc.
        """
        column_widths = self.column_widths()

        # Format the table headers. All of the items are centered
        total_length = sum(column_widths)

        # Prepare the talble for output
        table = ''
        # Add title, if present
        if isinstance(self.title, str):
            table += ''.join(('Table: ', self.title, '\n'))

        # Add top bar
        table += '\n' + '-' * total_length + '\n'

        # Add headers
        table += print_lines(self.column_titles, column_widths)

        table += '\n'

        # Add header bottom bars
        table += ''. join(['-' * (width)
                           if count == 0 else ' ' + '-' * (width - 1)
                           for count, width in enumerate(column_widths)])
        table += '\n'

        # Add rows
        for row in self.rows:
            table += print_lines(row, column_widths)
            table += ('\n\n' if self.multiline else '\n')

        # Add bottom bar
        table += '-' * total_length + '\n\n'

        return table
