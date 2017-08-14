"""
Test the FormattedStr class
"""

import unittest

from mollib.utils import FormattedStr


class TestFormattedStr(unittest.TestCase):

    def test_formatted_str(self):

        # Turn off the colored terminal
        FormattedStr.color_term = False
        self.assertEqual(FormattedStr('test'), 'test')
        self.assertEqual(FormattedStr('test', 'green'), 'test')
        self.assertEqual(len(FormattedStr('test', 'green')), 4)

        # Turn on the colored terminal
        FormattedStr.color_term = True
        self.assertEqual(FormattedStr('test'), 'test')
        self.assertNotEqual(FormattedStr('test', 'green'), 'test')

        # Disable strip_formatting_from_len
        FormattedStr.strip_formatting_from_len = False
        self.assertEqual(len(FormattedStr('test', 'green')), 13)

        # Enable strip_formatting_from_len
        FormattedStr.strip_formatting_from_len = True
        self.assertEqual(len(FormattedStr('test', 'green')), 4)

        # Try the justification functions
        self.assertEqual(len(FormattedStr('test', 'green').ljust(6)),
                         6)
        self.assertIsInstance(FormattedStr('test', 'green').ljust(6),
                              FormattedStr)
        self.assertEqual(len(FormattedStr('test', 'green').rjust(6)),
                         6)
        self.assertIsInstance(FormattedStr('test', 'green').rjust(6),
                              FormattedStr)
        self.assertEqual(len(FormattedStr('test', 'green').center(6)),
                         6)
        self.assertIsInstance(FormattedStr('test', 'green').center(6),
                              FormattedStr)