"""MolLib command line interface
"""
# Author: Justin L Lorieau
# Copyright 2016

import argparse
import mollib

def main():
    parser = argparse.ArgumentParser(prog='mollib',
                                     description='A molecular processor')
    parser._optionals.title = 'options'

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {}'.format(mollib.__version__))
    parser.add_argument('-i', '--in', action='append', nargs='+',
                        type=str, required=True, metavar='identifier',
                        help=('(required) The filename(s) or PDB identifier(s)'
                              ' of the structure(s)'))
    parser.add_argument('-o', '--out', nargs=1,
                        type=str, required=False, metavar='output_file',
                        help='The output filename of the structure')
    parser.add_argument('-c', '--config', nargs=1,
                        type=str, required=False, metavar='config_file',
                        help='The configuration file to use')
    parser.add_argument('--hydrogenate', action='store_true',
                        help='Strip the hydrogens and re-add hydrogens')
    parser.add_argument('--hbonds', required=False,
                        choices=['all', 'amide', 'aliphatic'],
                        help='Determine the hydrogen bonds.3')

    args = parser.parse_args()

if __name__ == "__main__":
    main()