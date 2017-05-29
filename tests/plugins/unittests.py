"""Unittests for the plugin module.
"""

import unittest
import argparse

from mollib.plugins import PluginManager


class TestMolLib(unittest.TestCase):

    def test_basic_plugin_subclass(self):
        """Tests basic subclasses of plugins. These tests include the Process
        and measure plugins, which are part of the core."""

        # Test the created subparser
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers(title='commands')

        # Setup the plugin manager
        plugin_manager = PluginManager(parser=parser, subparser=subparsers)
        parser = plugin_manager.process_parsers()

        # Check that the plugins were properly loaded by checking the 'process'
        # and 'measure' plugins
        self.assertIn('process',
                      [p.command for p in plugin_manager.plugins()])
        self.assertIn('measure',
                      [p.command for p in plugin_manager.plugins()])
