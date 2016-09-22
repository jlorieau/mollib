import unittest
import os

from mollib.statistics import Statistics, RamachandranStatistics


class TestStatistics(unittest.TestCase):
    "Test the Statistics class."

    def test_load(self):
        "Test the proper loading of the input files."
        stats = Statistics()

        # Make sure the identifiers are properly loaded.
        self.assertGreater(len(stats.identifiers), 1)
        self.assertIsInstance(stats.identifiers.pop(), str)

        # Make sure the output filenames are correct and correspond to real
        # paths
        filename =  stats.get_measurement_filename()
        path, filename = os.path.split(filename)
        self.assertTrue(os.path.isdir(path))


