"""
MolLib default settings
"""
# Author: Justin L Lorieau
# Copyright: 2016

# Multiple implementations are possible here. I decided to put the settings
# in the module itself because the settings are then documented well by Sphinx.
# An alternative implementation that returns a singleton dict would not
# document well. Returning a settings object would also be harder to access
# through the api.
#
# In the implementation selected, the integrity of the module is kept by the
# import_settings function, which makes sure that the config setting
# matches that within the module. Also, functions cannot be replaced.

import ast
import logging
from collections import OrderedDict


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

# Create the _settings_modules dict and register the core settings module
_setting_modules = OrderedDict({'settings': globals()})


def register_settings(module):
    """Registers a settings module.

    The name of the module is inferred from the module.__name__ field.
    """
    global _setting_modules
    module_name = module.__name__.strip('mollib.')

    if module_name in _setting_modules:
        msg = "The settings for '{}' are already registered."
        logging.error(msg.format(module_name))
    _setting_modules[module_name] = module


def list_global_settings():
    """Lists the settings modules currently installed.

    Returns
    -------
    list of str
        A list of the settings strings currently installed. The names returned
        are those that can be directly edited in configuration file sections.
    """
    global _setting_modules
    return _setting_modules.keys()


def import_config(config):
    """Import all of the settings blocks from the configparser object.

    Parameters
    ----------
    config: ``configparser.ConfigParser``
        The ConfigParser object with all of the custom settings.
    """
    # # Import this module's settings
    # import_settings(config, 'settings', locals())

    # Import other module settings in the core
    global _setting_modules

    for setting_name, setting_module in _setting_modules.items():
        import_settings(config, setting_name, setting_module)


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


    .. note:: The config parser only stores case-insensitive parameter names.
              Consequently, setting parameters should not have the same
              lowercase names.
    """
    # TODO: add string type validation

    err_msg = ("Parameter '{param}' in section '{sect}' could "
               "not be read.\n Expected a '{type_expected}'.")

    # If the settings_dict is a module instead of a dict, get the module's
    # dict instead
    if (not isinstance(settings_dict, dict) and
        hasattr(settings_dict, '__dict__')):
        settings_dict = settings_dict.__dict__

    # Find all of the parameters in the settings_dict and replace them with
    # values found in the config, provided the parameters are not callable
    # functions.
    for parameter, value in settings_dict.items():
        if config.has_option(section, parameter) and not callable(value):

            # Check the parameter value's type and convert accordingly
            # These must match. If they don't an error is raise, and we
            # proceed to the next value
            if type(value) == int:
                try:
                    config_value = config.getint(section, parameter)
                except ValueError:
                    logging.error(err_msg.format(param=parameter, sect=section,
                                                 type_expected='int'))
                    continue

            elif type(value) == float:
                try:
                    config_value = config.getfloat(section, parameter)
                except ValueError:
                    logging.error(err_msg.format(param=parameter, sect=section,
                                                 type_expected='float'))
                    continue

            elif type(value) == bool:
                try:
                    config_value = config.getboolean(section, parameter)
                except ValueError:
                    logging.error(err_msg.format(param=parameter, sect=section,
                                                 type_expected='bool'))
                    continue

            else:
                # The executable (callable) parameters cannot be replaced.
                # See the top of the if statement.
                config_value = config.get(section, parameter)  # str
                config_value = ast.literal_eval(config_value)

                # The config_value has to match the same type as the value it
                # seeks to replace. Moreover, it must be a base type, like a
                # string, dict, set, list or tuple.
                if (type(value) != type(config_value) or
                    not any(isinstance(config_value, i)
                            for i in [str, dict, set, list, tuple])):
                    logging.error(err_msg.format(param=parameter, sect=section,
                                                 type_expected=type(value)))
                    continue

            # Set the value
            settings_dict[parameter] = config_value

            # Debug message
            msg = "Config section '{}' imported '{}'"
            logging.debug(msg.format(section, parameter))
