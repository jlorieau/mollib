"""MolLib command line interface
"""
# Author: Justin L Lorieau
# Copyright 2016

import argparse

import mollib
from mollib.plugins import Plugin


def main():
    # Load the plugins
    plugins = Plugin.plugin_instances()

    # Load the argument parser and subparsers
    parser = argparse.ArgumentParser(prog='mollib',
                                     description='A molecular processor')
    subparsers = parser.add_subparsers(title='commands', dest='command',
                                       metavar='')
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' + mollib.__version__))

    # Process the plugin subparsers
    for plugin in plugins:
        plugin.options(subparsers)

    # Parse the commands
    args = parser.parse_args()

    # Prepare and preprocess the structure
    molecules = [mollib.Molecule(identifier) for identifier in args.i[0]]
    print(args)

    # Find the relevant plugins to execute
    active_plugins = [plugin for plugin in plugins if plugin.selected(args)]

    for molecule in molecules:
        for plugin in active_plugins:
            plugin.process(molecule)


if __name__ == "__main__":
    main()