"""
The settings manager functions to register and load settings.

Settings are stored in 'settings.py' files for each submodule or plugin.
Accesses to values in the settings should be made directly, rather than copied,
in the code so that updated settings (from config files) can be incorporated
properly.
"""

import ast
import os
import logging
from collections import OrderedDict

from . import settings

# Create the _settings_modules dict and register the core settings module
_settings_modules = OrderedDict({'settings': settings})


def register_settings(module):
    """Registers a settings module.

    The name of the module is inferred from the module.__name__ field.
    """
    global _settings_modules
    module_name = module.__name__.strip('mollib.')

    if module_name in _settings_modules:
        msg = "The settings for '{}' are already registered."
        logging.error(msg.format(module_name))
    _settings_modules[module_name] = module


def list_global_settings():
    """Lists the settings modules currently installed.

    Returns
    -------
    list of str
        A list of the settings strings currently installed. The names returned
        are those that can be directly edited in configuration file sections.
    """
    global _settings_modules
    return _settings_modules.keys()


def load_settings(config=None):
    """Import all of the settings blocks from the configparser object.

    Parameters
    ----------
    config: ``configparser.ConfigParser`` or None
        The ConfigParser object with all of the custom settings.
        Or None.
    """
    # Import other module settings in the core
    global _settings_modules

    for setting_name, setting_module in _settings_modules.items():
        import_settings(config, setting_name, setting_module)


def import_settings(config, section, settings_dict):
    """Import a settings block from the configparser object.

    This function works in concert with load_settings.

    Parameters
    ----------
    config: ``configparser.ConfigParser`` or None
        The ConfigParser object with all of the custom settings.
        Or None.
    section: str
        The name of the settings section to import.
    settings_dict: dict
        The dict for the settings module.


    .. note:: The config parser only stores case-insensitive parameter names.
              Consequently, setting parameters should not have the same
              lowercase names.

    .. note:: Any settings with the string 'path' in it will have the root
              path appended to it.
    """
    # Get the root path for the application
    root_path = os.path.abspath(os.path.dirname(__file__))

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

        # Load the configuration, if is is specified)
        if (config and config.has_option(section, parameter) and
           not callable(value)):

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

        # Any parameter with 'path' in it shall of the root path added to it
        if 'path' in parameter:
            settings_dict[parameter] = os.path.join(root_path, '..', value)