"""
Utility functions for dealing with numbers.
"""
from math import floor, log10


def center(number, rjust=4):
    """Center a number string about the decimal."""
    number = str(number)
    split = number.split('.')
    return (''.join((split[0].rjust(4), '.', split[1]))
            if len(split) > 1 else split[0])


def round_sig(x, sig=2):
    """Round a number to the specified number of significant figures."""
    return round(x, sig-int(floor(log10(abs(x)))) - 1)
