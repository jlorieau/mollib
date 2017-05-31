Overview
========

The command line interface includes all of mollib's functions for processing
and analyzing molecules and molecular data. Mollib can be accessed either
through the ``mollib`` or ``ml`` commands.

.. include:: cmds/ml_help.rst


``-h`` / ``--help``
    Get basic help information on usage of the command line interface, or one
    one of the mollib commands.

``-d`` / ``--debug``
    Display debug messages to the terminal. This is usually only useful to
    developers or in helping to debug issues.

``-s`` / ``--suppress``
    Display only errors and critical errors.

``-v`` / ``--verbose``
    Display informational messages to the terminal. By default, only warning
    and error messages are displayed.

``--version``
    Display the version number of the installed mollib.

``--list-plugins``
    Display a list of installed and enabled plugins.

    .. include:: cmds/ml_list_plugins.rst

``--list-settings``
    Display a list of the settings sections that will be interpreted from
    configuration files.

    .. include:: output/ml_list-settings.rst

.. _configuration-files:

Configuration Files
*******************

Mollib will check for a ``.mollibrc`` configuration file and a configuration
file passed with the ``-c`` ``filename`` / ``--config`` ``filename``
commandline argument for configuration parameters.

Configuration files customize mollib's behavior from default values. If a
``.mollibrc`` file is present and a configuration file is specified in the
commandline, the commandline configuration file will take precedence for
parameters that are mentioned in both.

The following is an example configuration file. Note that only the parameters
that need to be changed can be specified in the file.

.. literalinclude:: ../../examples/setup.cfg
    :caption: setup.cfg

.. TODO: Add a recipes section
