"""Unittests for the plugin module.
"""

import unittest
import argparse

from mollib.plugins import Plugin


class TestMolLib(unittest.TestCase):

    def test_basic_plugin_subclass(self):
        "Tests basic subclasses of plugins"

        class Submethod(Plugin):
            pass

        plugin = Submethod()
        self.assertIn(plugin, Submethod.plugin_instances())

        # Test the created subparser
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers(title='commands')
        subsubparser = plugin.options(subparsers)

        # Test that the basic parameters were added
        self.assertIn('-h', subsubparser._option_string_actions)

        # Test the process_molecule method. The base class should raise
        # a NotImplementedError
        with self.assertRaises(NotImplementedError):
            plugin.process(molecule=None)

