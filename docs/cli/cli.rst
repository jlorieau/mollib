=====================
Commandline Interface
=====================
The command line interface includes all of the mollib functions for processing
molecules.

    .. include:: output/cli_help.html


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

        .. include:: output/cli_list-plugins.html

    ``--list-settings``
        Display a list of the settings sections that will be interpreted from
        configuration files.

        .. include:: output/cli_list-settings.html

********
Commands
********

Commands include various types of operations that can be conducted by mollib.
The command listing may contain additional commands, depending on whether
mollib plugins are installed.

.. toctree::
    :maxdepth: 1

    process
    measure
    hbonds

********
Settings
********

    The following are the general configuration file settings.

    .. literalinclude:: /../mollib/core/settings.py
        :caption: [settings]

    The ``utils.settings`` impact how the data are presented.

    .. literalinclude:: /../mollib/utils/settings.py
        :caption: [utils.settings]

    The ``hydrogens.settings`` impact how hydrogen atoms are added to a
    molecule. These impact the ``--hydrogenate``.

    .. literalinclude:: /../mollib/hydrogens/settings.py
        :caption: [hydrogens.settings]

    The ``hbonds.settings`` impact how hydrogen bonds are detected. These impact
    the ``hbonds`` command and the ``--rama`` option in the ``measure`` command.

    .. literalinclude:: /../mollib/hbonds/settings.py
        :caption: [hbonds.settings]

    The following options only impact the creation of mollib datasets through
    the ``build_data`` setup command.

    .. literalinclude:: /../mollib/statistics/settings.py
        :caption: [statistics.settings]
