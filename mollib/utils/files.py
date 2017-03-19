import logging
import os

from . import settings


def check_file(filename, critical=True):
    """Check whether the filename exists.

    Paramaters
    ----------
    critical: bool (optional)
        If the existence of the file or path is critical to proceed, and the
        file doesn't exist, end the program execution.

    Returns
    -------
    None
    """
    if not os.path.isfile(filename):
        msg = "Filename '{}' does not exist."
        logging.error(msg.format(filename))
        if critical:
            exit()
    return None


def write_txt_file(filename, string, overwrite=None):
    pass