"""MolLib command line interface
"""
# Author: Justin L Lorieau
# Copyright 2016

import argparse
import logging
import sys
import os

import mollib
from mollib.utils import FormattedStr
from mollib.plugins import PluginManager
from mollib.core import list_global_settings, load_settings
import mollib.utils.settings


try:
    import configparser
except ImportError:
    import ConfigParser as configparser


def list_plugins(plugin_manager):
    "Prints a list of the installed plugins."
    print('Installed plugins:')
    for plugin in plugin_manager.plugins():
        msg = '\t{:<15} '.format(plugin.name)
        enabled = (FormattedStr('Enabled', 'green') if plugin.enabled else
                   FormattedStr('Not Enabled', 'red'))
        print(msg + enabled)


def list_settings():
    "Prints a list of the installed setting sections."
    print('Installed settings sections:')
    for section in list_global_settings():
        msg = '\t[{}]'.format(section)
        print(msg)


def main():
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
                        version=('%(prog)s ' + mollib.__version__.__version__),
                        help='Show the program version')

    # Load the plugins
    plugin_manager = PluginManager(parser=parser, subparser=subparsers)
    parser = plugin_manager.process_parsers()

    # process the --list-settings and --list_plugins options
    if '--list-plugins' in sys.argv:
        list_plugins(plugin_manager)
        exit()
    if '--list-settings' in sys.argv:
        list_settings()
        exit()

    # Parse the commands
    args = parser.parse_args()

    # Set special flags that need to be set before processing molecules
    if args.save:
        mollib.utils.settings.save_fetched_files_locally = True

    # Setup the logger
    fmt = '{}: %(levelname)-8s %(message)s'.format('mollib')
    logging.basicConfig(format=fmt, level=args.loglevel)
    logging.debug(args)

    # Read in the configuration file(s)
    config_files = [os.path.expanduser('~/.mollibrc'), ]
    if args.config:
        config_files.append(args.config[0])
    config = configparser.ConfigParser()
    config.read(config_files)
    load_settings(config)

    # Load the molecules
    if args.models:
        molecules = []
        mr = mollib.MoleculeReader()
        model_ids = args.models
        for identifier in args.i[0]:
            molecules += mr.read(identifier, model_ids=model_ids)
    else:
        molecules = [mollib.Molecule(identifier) for identifier in args.i[0]]

    # Find the relevant plugins to execute
    active_plugins = plugin_manager.plugins()

    # Pre-process the molecules
    for plugin in active_plugins:
        logging.debug('Preprocessing:{}'.format(plugin))
        plugin.preprocess(molecules, args)

    # Process the molecules
    for plugin in active_plugins:
        logging.debug('Processing: {}'.format(plugin))
        plugin.process(molecules, args)

    # Post-process the molecules
    for plugin in active_plugins:
        logging.debug('Post-rocessing: {}'.format(plugin))
        plugin.postprocess(molecules, args)

if __name__ == "__main__":
    main()