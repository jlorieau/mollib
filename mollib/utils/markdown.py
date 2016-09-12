"""
Utilities for rendering information in Markdown
"""

class MDTable(object):
    """Renders a table in Markdown.

    Attributes
    ----------
    title: str
        If specified, this attribute will be used as the table's title.
    column_titles: iterable
        An iterable (list) of strings for the column header titles
    rows: iterable
        An iterable (list) of rows. Each row is an iterable (list) of strings
        that match the number of columns in column_titles
    multiline: bool
        If True, this table will render as a multiline table
    """

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
        self.multiline = False
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
        "Return the list of column widths for the current data."
        # Find the largest item for each column
        column_widths = [len(c) for c in self.column_titles]

        for col_index in range(len(column_widths)):
            max_width = column_widths[col_index]

            for row in self.rows:
                if len(row[col_index]) > max_width:
                    max_width = len(row[col_index])
            column_widths[col_index] = max_width
        return column_widths

    def table_width(self):
        "Return the width of the table."
        column_widths = self.column_widths()
        return sum([width + 3 for width in column_widths]) - 1

    def content(self):
        """Renders a string for a markdown table.

        Note that this renders multiline tables in pandoc.
        """
        column_widths = self.column_widths()

        # Format the table headers. All of the items are centered
        total_length = self.table_width()


        # Prepare the talble for output
        table = ''
        # Add title, if present
        if isinstance(self.title, str):
            table += ''.join(('Table: ', self.title, '\n'))

        # Add top bar
        table += '\n' + '-' * total_length + '\n'

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
            table += ('\n\n' if self.multiline else '\n')

        # Add bottom bar
        table += '-' * total_length + '\n\n'

        return table
