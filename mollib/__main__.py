"""MolLib command line interface
"""
# Author: Justin L Lorieau
# Copyright 2016

import argparse
import mollib

def main():
    parser = argparse.ArgumentParser(prog='mollib',
                                     description='A molecular processor')
    # parser._optionals.title = 'arguments'

    subparsers = parser.add_subparsers(title='commands')

    # Commands and subparsers
    # Measurement subparser
    mparser = subparsers.add_parser('measure',
                                    help='Measure geometries of molecules')
    mparser.add_argument('-md', '--measure-distance', nargs=2, required=False,
                        metavar='atom',
                        help=("Measured distance between atoms.'"
                              "ex: N31-N C32-CA"))

    hbparser = subparsers.add_parser('hbond',
                                     help='Find hydrogen-bonds of molecules')
    hbparser.add_argument('-t', '--type', required=True,
                          choices=['all', 'amide', 'aliphatic'],
                          help='Determine the hydrogen bonds.')

    # Standard commands
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {}'.format(mollib.__version__))
    parser.add_argument('id', action='append', nargs='+',
                        type=str, metavar='identifier/filename',
                        help=('The filename(s) or PDB identifier(s)'
                              ' of the structure(s)'))
    parser.add_argument('-o', '--out', nargs=1, required=False,
                        type=str, metavar='filename',
                        help='The output filename of the structure')
    parser.add_argument('-c', '--config', nargs=1, required=False,
                        type=str, metavar='filename',
                        help='The configuration filename to use')

    # Processing options
    parser.add_argument('--hydrogenate', action='store_true',
                        help='Strip the hydrogens and re-add hydrogens')

    args = parser.parse_args()

if __name__ == "__main__":
    main()