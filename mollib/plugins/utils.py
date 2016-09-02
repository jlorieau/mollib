import argparse

def check_number_arguments(min_nargs, max_nargs, parameters_name=None):
    """Return a subclass of CheckNumberArguments with the min_nargs and
    max_nargs set.

    This function is used for the 'action' parameter of the argparse
    add_argument function to assert that the number of arguments pass fits
    within a range.

    Parameters
    ----------
    min_nargs: int
        The minimum number of arguments required.
    max_nargs: int
        The maximum number of arguments required.
    parameters_name: str (optional)
        If specified, this string will be used in error messages.
        ex: 'atoms' would display "At least 4 atoms should be specified.

    Returns
    -------
    class
        A subclass of CheckNumberArguments class
    """
    class CheckNumberArgumentsSubclass(CheckNumberArguments):
        pass
    CheckNumberArgumentsSubclass.parameters_name = parameters_name
    CheckNumberArgumentsSubclass.min_nargs = min_nargs
    CheckNumberArgumentsSubclass.max_nargs = max_nargs
    return CheckNumberArgumentsSubclass


class CheckNumberArguments(argparse.Action):
    "An argparse.Action class to test that the number of parameters was passed."

    parameters_name = None

    def __call__(self, parser, args, values, option_string=None):
        if len(values) < self.min_nargs:
            msg = '{}: At least {} {} should be specified'
            raise argparse.ArgumentError(None,
                                         msg.format(option_string,
                                                    self.min_nargs,
                                                    self.parameters_name or
                                                    'parameters'))
        if len(values) > self.max_nargs:
            msg = '{}: At most {} {} should be specified'
            raise argparse.ArgumentError(None,
                                         msg.format(option_string,
                                                    self.max_nargs,
                                                    self.parameters_name or
                                                    'parameters'))
        setattr(args, self.dest, values)
