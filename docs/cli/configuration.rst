Configuration Options
=====================

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
