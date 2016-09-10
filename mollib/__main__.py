"""MolLib command line interface
"""
# Author: Justin L Lorieau
# Copyright 2016

import argparse
import logging
import sys
import os

import mollib
from mollib.plugins import Plugin
from mollib.core.settings import list_global_settings, import_config

try:
    import configparser
except ImportError:
    import ConfigParser as configparser


def list_plugins():
    "Prints a list of the installed plugins."
    print('Installed plugins:')
    for plugin in Plugin.plugin_instances():
        msg = '\t{:<15} '.format(plugin.name)
        enabled = ('(\033[92mEnabled\033[0m)' if plugin.enabled
                   else '(\033[91mNot Enabled\033[0m)')
        print(msg + enabled)


def list_settings():
    "Prints a list of the installed setting sections."
    print('Installed settings sections:')
    for section in list_global_settings():
        msg = '\t[{}]'.format(section)
        print(msg)


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
                       help="Suppress all messages, except critical", )
    group.add_argument('-v', '--verbose',
                        action="store_const", dest="loglevel",
                        const=logging.INFO,
                        help="Print extra information")

    # Version information and other installation information
    parser.add_argument('--list-plugins',
                       action='store_true',
                       help='List the installed plugins')
    parser.add_argument('--list-settings',
                       action='store_true',
                       help='List the available sections for settings')
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' + mollib.__version__),
                        help='Show the program version')

    #TODO: Add argument to just download pdb file ('-g --get')

    # Process the plugin subparsers
    for plugin in plugins:
        plugin.options(subparsers)

    # process the --list-settings and --list_plugins options
    if '--list-plugins' in sys.argv:
        list_plugins()
        exit()
    if '--list-settings' in sys.argv:
        list_settings()
        exit()

    # Parse the commands
    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)
    logging.debug(args)

    # Read in the configuration file(s)
    config_files = [os.path.expanduser('~/.mollibrc'), ]
    if args.config:
        config_files.append(args.config)
    config = configparser.ConfigParser()
    config.read(config_files)
    import_config(config)

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