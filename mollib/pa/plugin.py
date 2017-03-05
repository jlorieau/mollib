# -*- coding: utf-8 -*-
"""
The plugin for the hbond submodule.
"""
from mollib.plugins import Plugin


class PA(Plugin):
    """The core plugin to offer the 'pa' command."""

    enabled = True
    order = 110
    create_command_subparser = True

    def help(self):
        return "Fit partial alignment RDC and RACS data using SVD"

    def process_parser(self):
        "Process the parser for the 'pa' command."
        p = self.command_subparsers['pa']

        p.add_argument('-a', '--alignment',
                       action='store_true',
                       required=True,
                       help="Alignment file with RDC and RACS data")

        p.add_argument('-p', '--pred',
                       action='store_true',
                       help="Report predicted RDCs and RACS")

