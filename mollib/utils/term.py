"""
Utility functions for terminals.

Parameters
----------
:terminal: bool
    Whether the program is being run in the terminal.
:rows: int
    The number of rows available in the terminal.
:columns: int
    The number of columns available in the terminal.
"""

import os
import sys

if sys.stdout.isatty():
    terminal = True
    rows, columns = os.popen('stty size', 'r').read().split()
else:
    terminal = False
    rows, columns = None, None
