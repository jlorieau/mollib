Hbonds Command
==============
The ``hbonds`` command detects and reports hydrogen bonds in molecules.

    .. include:: output/cli_hbonds_help.html

    **Examples**

    .. include:: output/cli_hbonds_i_2KXA.html

Arguments
---------

    ``--aliphatic``
        Include aliphatic hydrogen bonds in the results. The acceptor-donor
        atom cutoff distances are elongated to 3.0 A, and carbon atoms are
        allowed in hydrogen bond donor dipoles.

    ``--detailed``
        Present a detailed report on the geometries of hydrogen bonds.

Configuration File Settings
---------------------------
The following are the configuration file settings. Do not include the comments
in the file.

    .. literalinclude:: /../mollib/hbonds/settings.py
        :caption: [hbonds.settings]
