"""
MolLib default settings
"""
# Author: Justin L Lorieau
# Copyright: 2016

import ast
import logging


# Ramachandran parameters

dihedral_helix_angles = (-60, -45)
"The default optimal Ramachandran angles for a helix (in degrees)."

dihedral_helix_angles_threshold = (30, 30)
"""The default Ramanchandran angle threshold (+/- threshold) to be considered a
helix (in degrees)"""

default_pH = 7.0
"""The default pH of new molecules."""

pKs = {'ASP': {'OD1-OD2': (-1.0, 3.5)},
       'GLU': {'OE1-OE2': (-1.0, 4.2)},
       'HIS': {'ND1-NE2': (6.6, 14.0)},
       'CYS': {'SG': (6.8,)},
       'TYR': {'OH': (10.3,)},
       'LYS': {'NZ': (10.5, 14.0, 14.0)},
       'last': {'O-OXT': (-1.0, 3.3)},
       'first': {'N': (7.7, 14.0, 14.0),}
     }
"""The default pKs of ionizable amino-acids. [Ref]_

Some amino-acids have degenerate ionizeable atoms; these are listed and
separated by '-' characters. The different pKs for each ionization is
listed in the subsequent items in the tuple.


  .. [Ref] G. R. Grimsley, J. M. Scholtz, C. N. Pace, Protein Sci. 18, 247-51
           (2009).
"""


def import_config(config):
    """Import all of the settings blocks from the configparser object.

    Parameters
    ----------
    config: ``configparser.ConfigParser``
        The ConfigParser object with all of the custom settings.
    """
    # Import this module's settings
    import_settings(config, 'settings', locals())

    # Import other module settings
    setting_sections = [s for s in config.sections()
                        if s.startswith('settings.')]
    # for section in setting_sections:
    #     try

def import_settings(config, section, settings_dict):
    """Import a settings block from the configparser object.

    This function works in concert with import_config.

    Parameters
    ----------
    config: ``configparser.ConfigParser``
        The ConfigParser object with all of the custom settings.
    section: str
        The name of the settings section to import.
    settings_dict: dict
        The dict for the settings module.
    """
    err_msg = ("Parameter '{param}' in section '{sect}' could "
               "not be read.\n Expected a '{type_expected}'.")

    # Find all of the parameters in the settings_dict and replace them with
    # values found in the config, provided the parameters are not callable
    # functions.
    for parameter, value in settings_dict.items():
        if config.has_option(section, parameter) and not callable(parameter):

            # Check the parameter value's type and convert accordingly
            # These must match. If they don't an error is raise, and we
            # proceed to the next value
            if type(value) == int:
                try:
                    config_value = config.getint('settings', parameter)
                except ValueError:
                    logging.error(err_msg.format(param=parameter, sect=section,
                                                 type_expected='int'))
                    continue

            elif type(value) == float:
                try:
                    config_value = config.getfloat('settings', parameter)
                except ValueError:
                    logging.error(err_msg.format(param=parameter, sect=section,
                                                 type_expected='float'))
                    continue

            elif type(value) == bool:
                try:
                    config_value = config.getboolean('settings', parameter)
                except ValueError:
                    logging.error(err_msg.format(param=parameter, sect=section,
                                                 type_expected='bool'))
                    continue

            else:
                config_value = config.get('settings', parameter)  # str
                config_value = ast.literal_eval(config_value)
                if type(value) != type(config_value):
                    logging.error(err_msg.format(param=parameter, sect=section,
                                                 type_expected=type(value)))
                    continue

            # Set the value
            settings_dict[parameter] = config_value

            # Debug message
            msg = "Config section '{}' imported '{}'"
            logging.debug(msg.format(section, parameter))
