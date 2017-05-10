"""
Utilities for rendering information in Markdown
"""
from math import ceil

from . import settings
from . import term
from .formatted_str import FormattedStr, wrap


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
    >>> t
    'My             first and ever multiline text '
    >>> t = print_lines(['My', 'first and ever', 'multiline text'], \
                        [10, 10, 10])
    >>> t.splitlines()
    ['My        first and multiline ', '          ever      text      ']
    """
    assert(len(text_items) == len(widths))

    # wrap the text_items
    text_items = [wrap(t, w) for t,w in zip(text_items, widths)]
    no_lines = [len(t) for t in text_items]
    max_no_lines = max(no_lines)
    text = ''
    for line_no in range(max_no_lines):
        for text_item, width in zip(text_items, widths):
            text += (text_item[line_no].ljust(width)
                     if line_no < len(text_item) else
                     ''.ljust(width))

        # Print a new line on all lines except the last one
        if line_no != max_no_lines - 1 :
            text += '\n'
    return text


class MDTable(object):
    """Renders a table in Markdown.

    Attributes
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
    max_width = None

    def __init__(self, *args):
        """Constructor that initializes the column headers.

        The current implementation only uses single-line rows and headers.
        (No column columnwidth maximum)

        Parameters
        ----------
        *args: iterable
           A iterable (list) of strings for the column header titles
        """
        self.column_titles = [FormattedStr(i, 'bold') for i in args]
        self.rows = []
        self.title = None

        if term.terminal and term.columns is not None:
            self.max_width = int(term.columns)
        else:
            self.max_width = settings.default_max_width

    @property
    def empty_headers(self):
        """True if the headers are all empty (i.e. equal to '')."""
        stripped_headers = [i.stripped_str() if hasattr(i, 'stripped_str')
                            else i
                            for i in self.column_titles]
        return all([i == '' for i in stripped_headers])

    def num_cols(self):
        "Return the number of columns in the table."
        return len(self.column_titles)

    def add_row(self, *args):
        """Add a row of values to the table"""
        # The number of items in the row must match the number of columns
        assert(len(args) == len(self.column_titles))

        # Add the rows and convert them to strings
        self.rows.append([str(i) if not isinstance(i, str) else i
                         for i in args])

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

                # Resort the largest_column_widths by index number, and copy
                # over to the column widths
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

        # Prepare the table for output
        table = ''
        # Add title, if present
        if isinstance(self.title, str):
            table += FormattedStr('Table: ', 'bold') + self.title + '\n'

        # Add top bar. Only needed for multiline tables that have headers
        if self.multiline and not self.empty_headers:
            table += '\n' + '-' * total_length + '\n'
        else:
            table += '\n'

        # Add headers
        if not self.empty_headers:
            table += print_lines(self.column_titles, column_widths)
            table += '\n'

        # Add header bottom bars
        broken_lines = ''. join(['-' * (width-1)
                           if count == 0 else ' ' + '-' * (width-1)
                           for count, width in enumerate(column_widths)])
        table += broken_lines
        table += '-\n'

        # Add rows
        for row in self.rows:
            table += print_lines(row, column_widths)
            table += ('\n\n' if self.multiline else '\n')

        # Add bottom bar. Only needed for multiline tables or if the column
        # headers are all empty
        if self.multiline:
            table += '-' * total_length + '\n'
        elif self.empty_headers:
            table += broken_lines
            table += '-\n'
        return table


def dict_table(dictionary, sort_key=None):
    """Renders a Markdown table for a dictionary.

    Parameters
    ----------
    dictionary: dict
        The dictionary to prepare in the table.
    sort_key: function
        If set, the dict keys will be sorted by the given sort function.

    Returns
    -------
    table: :obj:`mollib.utils.MDTable`
        The table generated from the dictionary.
    """
    # Determine the number of columns. One column for the keys, and one for
    # each item in an iterable of the dict values
    lengths = [len(v) if hasattr(v, '__len__') and not isinstance(v, str) else 1
               for v in dictionary.values()]
    no_cols = max(lengths) + 1

    # Create the table with empty headers
    table = MDTable(*['' for i in range(no_cols)])

    # Get the dict keys to start populating the rows
    if sort_key is not None:
        keys = sorted(dictionary.keys(), key=sort_key)
    else:
        keys = dictionary.keys()

    # Populate the rows
    for key in keys:
        values = dictionary[key]

        # Fill the list of values to match the number of columns
        if hasattr(values, '__len__') and not isinstance(values, str):
            # If the values are dicts, then use these to list the key/value
            # pairs. Otherwise, just list the values themselves

            if isinstance(values, dict):
                values = ["{}: {}".format(k, v)
                          for k,v in values.items()]
            values = list(values) + [''] * (no_cols - len(values) - 1)
        else:
            values = [values, ] + [''] * (no_cols - 2)

        table.add_row(key, *values)

    return table




