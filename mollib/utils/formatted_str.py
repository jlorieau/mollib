"""
Classes and functions to format strings to various output formats.
"""
import sys
import os
import re
import textwrap

from . import term


re_ansi = re.compile(r'(\x9B|\x1B\[)[0-?]*[ -\/]*[@-~]')


def wrap(string, *args, **kwargs):
    "A text wrapper than properly works with str and FormattedStr objects."
    if isinstance(string, FormattedStr) and string.formatters is not None:
        # First strip the formatting information
        stripped_string = string.stripped_str()

        # Then wrap the string and rewrap all of the items in the formatted
        # tags.
        stripped_strings = textwrap.wrap(stripped_string, *args, **kwargs)
        stripped_strings = [FormattedStr(s, *string.formatters)
                            for s in stripped_strings]
        return stripped_strings
    else:
        return textwrap.wrap(string, *args, **kwargs)


class FormattedStr(str):
    """A formatted string.

    Attributes
    ----------
    color_term: bool, optional
        (class attribute) The string output is for a terminal that supports
        color.
    bw_term: bool, optional
        (class attribute) The output string is for a terminal that only support
        black and white, including PIPE output.
    html: bool, optional
        (class attribute) The output should be formatted html.
    convert_greek: bool, optional
        (class attribute) Textual greek (i.e.: 'alpha', 'beta', 'Phi') should
        be converted to their unicode Greek characters.
    strip_formatting_from_len: bool, optional
        (class attribute) If true, formatting text characters will not
        contribute to the length reported by __len__.
    formatters: tuple of str
        The formatters to use in formatting this string.


    .. note:: The FormattedStr prints colored ansi output if the terminal
              is capable of it or if the 'FORCE_COLOR' environment variable is
              set.

    .. note:: Output to PIPE strips colors.

    .. note:: FormattedStr turn into regular strings when joined with regular
              strings.


    Examples
    --------
    >>> t = FormattedStr('test')
    >>> print(t)
    test
    >>> len(t)
    4
    """

    color_term = False
    convert_greek = True
    strip_formatting_from_len = True

    formatters = None

    def __new__(cls, string, *args):
        """The constructor

        Parameters
        ----------
        args: tuple of str
            The first item is the string to be formatted. The subsequent items
            are formatters (str) on how to format this string. Available
            formatting options include: 'red', 'green', 'yellow', 'blue',
            'magenta', 'cyan', 'bold', 'underline'.
        """
        if cls.color_term:
            string = cls.color_term_formatter(string, *args)
        fstring =  super(FormattedStr, cls).__new__(cls, string)
        fstring.formatters = args
        return fstring

    def __len__(self):
        if self.strip_formatting_from_len:
            return self._stripped_len()
        return str.__len__(self)

    def _stripped_len(self):
        "Return the length of the stripped string."
        stripped_str = str(self)
        stripped_str = re_ansi.sub('', stripped_str)
        return str.__len__(stripped_str)

    def ljust(self, width, *args, **kwargs):
        diff = str.__len__(self) - self._stripped_len()
        new_width = width + diff
        return FormattedStr(str.ljust(self, new_width, *args, **kwargs))

    def rjust(self, width, *args, **kwargs):
        diff = str.__len__(self) - self._stripped_len()
        new_width = width + diff
        return FormattedStr(str.rjust(self, new_width, *args, **kwargs))

    def center(self, width, *args, **kwargs):
        diff = str.__len__(self) - self._stripped_len()
        new_width = width + diff
        return FormattedStr(str.center(self, new_width, *args, **kwargs))

    @classmethod
    def color_term_formatter(cls, string, *formatters):
        """Formats the string for a colored terminal."""
        for formatter in formatters:
            if formatter == 'red':
                string = ''.join(('\033[91m', string, '\033[0m'))
            elif formatter == 'green':
                string = ''.join(('\033[92m', string, '\033[0m'))
            elif formatter == 'yellow':
                string = ''.join(('\033[33m', string, '\033[0m'))
            elif formatter == 'blue':
                string = ''.join(('\033[94m', string, '\033[0m'))
            elif formatter == 'magenta':
                string = ''.join(('\033[95m', string, '\033[0m'))
            elif formatter == 'cyan':
                string = ''.join(('\033[96m', string, '\033[0m'))

            if formatter == 'bold':
                string = ''.join(('\033[1m', string, '\033[22m'))
            if formatter == 'underline':
                string = ''.join(('\033[4m', string, '\033[24m'))

        return string

    def stripped_str(self):
        "Return the stripped string (i.e. without formatting codes)."
        stripped_str = str(self)
        stripped_str = re_ansi.sub('', stripped_str)
        return stripped_str


# Setup the terminal
if term.terminal or os.environ.get('FORCE_COLOR', False):
    FormattedStr.color_term = True
