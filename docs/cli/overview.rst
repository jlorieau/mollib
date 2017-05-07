Overview
========

The command line interface includes all of mollib's functions for processing
and analyzing molecules and molecular data.

    .. include:: output/cli_help.rst


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

        **Examples**

        .. include:: output/cli_list-plugins.rst

    ``--list-settings``
        Display a list of the settings sections that will be interpreted from
        configuration files.

        .. include:: output/cli_list-settings.rst

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

Configuration Options
*********************

    The following are the general configuration file settings.

    .. literalinclude:: /../mollib/core/settings.py
        :caption: [settings]

    The ``utils.settings`` impact how general functions behave and how the data
    are presented.

    .. literalinclude:: /../mollib/utils/settings.py
        :caption: [utils.settings]

    The ``hydrogens.settings`` impact how hydrogen atoms are added to a
    molecule. These impact the ``--hydrogenate`` option.

    .. literalinclude:: /../mollib/hydrogens/settings.py
        :caption: [hydrogens.settings]

    The ``hbonds.settings`` impact how hydrogen bonds are detected. These impact
    the ``hbonds`` command and the ``--rama`` option in the ``measure`` command.

    .. literalinclude:: /../mollib/hbonds/settings.py
        :caption: [hbonds.settings]

    The ``pa.settings`` impact how partial alignment (RDC and RACS) data are
    fit.

    .. literalinclude:: /../mollib/pa/settings.py
        :caption: [pa.settings]

    The following options only impact the creation of mollib datasets through
    the ``build_data`` setup command.

    .. literalinclude:: /../mollib/statistics/settings.py
        :caption: [statistics.settings]
