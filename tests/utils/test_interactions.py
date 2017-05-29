"""
Tests for the interaction label utility functions.
"""
import unittest
import re

from mollib.utils.interactions import _re_interaction_str


class TestInteractions(unittest.TestCase):

    def test_interaction_regex(self):
        """Test the regex matching of interaction labels."""

        re_interaction = re.compile(_re_interaction_str)

        # The following regexes should match
        good_labels = ('14N-H',
                       'A.14N-H',
                       'A.14N-13H',
                       'A.14N-H',
                       'A.41N-B.32H-N',
                       '31N-C-2',
                       '32N-C+1',
                       '13C',
                       'C.14N',
                       '13N-CA-H+1',
                       '13CA-HA#',
                       )

        # The following regexes should not match
        bad_labels = ('13',
                      '2.223',
                      'CA',
                      )

        # Check that all of the good labels match completely.
        for label in good_labels:
            match = re_interaction.search(label)

            self.assertIsNotNone(match)
            self.assertEqual(match.group(), label)

        # Check that all of the bad labels do not match.
        for label in bad_labels:
            match = re_interaction.search(label)

            self.assertIsNone(match)
