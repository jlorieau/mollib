"""MolLib command line interface
"""
# Author: Justin L Lorieau
# Copyright 2016

import argparse
import logging

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

    # Logging levels
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-d', '--debug',
                        action="store_const", dest="loglevel",
                        const=logging.DEBUG,
                        default=logging.WARNING,
                        help="Print debugging statements",)
    group.add_argument('-s', '--suppress',
                       action="store_const", dest="loglevel",
                       const=logging.CRITICAL,
                       default=logging.WARNING,
                       help="Suppress all messages, except critical.", )
    group.add_argument('-v', '--verbose',
                        action="store_const", dest="loglevel",
                        const=logging.INFO,
                        help="Print extra information")

    # Version information
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' + mollib.__version__),
                        help='Show the program version')

    #TODO: Add argument to just download pdb file ('-g --get')

    # Process the plugin subparsers
    for plugin in plugins:
        plugin.options(subparsers)

    # Parse the commands
    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)
    logging.debug(args)
    logging.debug(plugins)

    # Prepare and preprocess the structure
    molecules = [mollib.Molecule(identifier) for identifier in args.i[0]]


    # Find the relevant plugins to execute
    active_plugins = [plugin for plugin in plugins if plugin.selected(args)]

    # Pre-process the molecules
    for molecule in molecules:
        for plugin in active_plugins:
            plugin.preprocess(molecule, args)

    # Process the molecules
    for molecule in molecules:
        for plugin in active_plugins:
            plugin.process(molecule, args)

    # Post-process the molecules
    for molecule in molecules:
        for plugin in active_plugins:
            plugin.postprocess(molecule, args)

if __name__ == "__main__":
    main()