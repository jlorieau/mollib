=====================
Commandline Interface
=====================
The command line interface includes all of the mollib functions for processing
molecules.

    .. literalinclude:: cli/output/cli_help.txt
        :language: shell-session

    ``-h`` / ``--help``
        Get basic help information on usage of the command line interface, or one
        one of the mollib commands.

    ``-d`` / ``--debug``
        Display debug messages to the terminal. This is usually only useful to
        developers or in helping to debug issues.

    ``-v`` / ``--verbose``
        Display informational messages to the terminal. By default, only warning
        and error messages are displayed.

    ``--version``
        Display the version number of the installed mollib.

********
Commands
********

Commands include various types of operations that can be conducted by mollib.
The command listing may contain additional commands, depending on whether
mollib plugins are installed.

.. toctree::
    :maxdepth: 1

    cli/process
    cli/measure
