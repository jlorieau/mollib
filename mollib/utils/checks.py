"""
Utilities to check conditions and print helpful messages. 

These are useful for reporting checks to the CLI, and they should only be used
in CLI functions and classes themselves.
"""

import logging
import inspect
import os


def check_file(filename, msg=None, critical=True):
    """Check whether the file at the given filename exists.

    Parameters
    ----------
    filename: str
        The name of the file to check.
    msg: str (optional)
        If specified, the following message will be logged if the check
        doesn't pass.
    critical: bool (optional)
        If the existence of the file or path is critical to proceed, and the
        file doesn't exist, end the program execution.

    Returns
    -------
    bool
        True, if the file exists, False if it doesn't.
    """
    if not os.path.isfile(filename):
        parent_stack = inspect.stack()[1]
        parent_func = parent_stack[3]
        parent_loc = parent_stack[1] + ":" + str(parent_stack[2])

        msg = ("Filename '{}' does not exist.".format(filename)
               if msg is None else msg)
        logging.error(parent_func + ": " + msg)
        logging.debug(parent_func + ": " + parent_loc)

        if critical:
            exit()
        return False
    return True

#: Add custom msg
def check_not_empty(msg=None, critical=True, **kwargs):
    """Check whether the kwargs objects are empty.

    Parameters
    ----------
    msg: str (optional)
        If specified, the following message will be logged if the check
        doesn't pass.
    critical: bool (optional)
        If True and the passed objects musn't be empty to proceed, then end
        program execution.
    kwargs: dict
        The items to be checked for emptiness.

    Returns
    -------
    bool
        True if *all* of the kwargs items are not empty,
        False if they are empty

    Examples
    --------
    >>> check_not_empty(critical=False, data={'a': 1})
    True
    >>> check_not_empty(critical=False, data={})
    False
    """
    # Parameters must be passed to the function
    assert len(kwargs) > 0

    not_empty = all([len(i) > 0 for i in kwargs.values()
                     if hasattr(i, '__len__')])

    if not_empty:
        return True
    else:
        parent_stack = inspect.stack()[1]
        parent_func = parent_stack[3]
        parent_loc = parent_stack[1] + ":" + str(parent_stack[2])

        for k,v in kwargs.items():
            if hasattr(v, '__len__') and len(v) == 0:
                msg = "No items in '{}'.".format(k) if msg is None else msg
                logging.error(parent_func + ": " + msg)
                logging.debug(parent_func + ": " + parent_loc)

        if critical:
            exit()
        return False
