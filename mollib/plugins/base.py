"""
The base Plugin class.
"""
import argparse
from abc import ABCMeta, abstractmethod


class PluginManager(object):
    """Manager for plugins.

    Manages subclass instances of :obj:`Plugin`.

    Attributes
    ----------
    parser: :obj:`argparser.ArgumentParser`
        The root parser.
    subparser: subparser
        The subparser for commands.
    """

    parser = None
    subparser = None

    def __init__(self, parser, subparser):
        self.parser = parser
        self.subparser = subparser

        # Create all of the subparsers for existing instances
        for instance in self.plugins():
            if instance.create_command_subparser:
                command = instance.command
                s = subparser.add_parser(command,
                                         help=instance.help())
                Plugin.command_subparsers[command] = s

    def process_parsers(self):
        """Process the parsers using all of the plugin instances and return the
        parser.

        The method will call all of the :meth:`Plugin.process_parser` methods for
        each plugin instance.

        Returns
        -------
        parser: :obj:`argparse.ArgumentParser`
            The processed parser.
        """
        # Process all of the parsers for sublcasses of this class.
        for instance in self.plugins():
            instance.process_parser()
        return self.parser

    def plugins(self):
        """Return the active plugin instances."""
        return Plugin._instances

class Plugin(object):
    """Abstract base class for plugins.

    Attributes
    ----------
    name: str, optional
        The name of the plugin.
    enabled: bool
        If true, the plugin is enabled.
    command: str, optional
        The shell command for this plugin.
    order: int
        The order to place the plugin in the command line options. (Higher is
        lower priority)
    create_command_subparser: bool, optional
        If True, a new command will be created for this plugin, base on the
        plugin subclasses's name (in lowercase).
    command_subparsers: dict
        A dict containing the command names (key) and the comand subparsers
        (value). Access this dict to add arguments to specific commands


    .. note:: Preprocessor and postprocessor plugins should have an order
              lower than 0 and should access the 'process' parent parser.
              Analysis methods should have an order between
              200-1000.
    """
    __metaclass__ = ABCMeta

    name = None
    enabled = True
    command = None
    order = 200

    create_command_subparser = True

    command_subparsers = {}

    def __new__(cls, *args, **kwargs):
        "Keep track of subclass instances"
        # Create the instance and store it in this Base classe's instance
        # list
        instance = object.__new__(cls, *args, **kwargs)
        if "_instances" not in Plugin.__dict__:
            Plugin._instances = []
        Plugin._instances.append(instance)
        return instance

    def __init__(self):
        self.name = self.__class__.__name__.lower()
        self.command = self.__class__.__name__.lower()

    @abstractmethod
    def process_parser(self):
        """Process the argument parser for this plugin.

        This function is used to add arguments to the argument parser or to
        the command subparsers.

        The subparser for a specific command is accessed by name from the
        self.command_subparsers. ex: self.command_subparsers['hbonds']

        This function is called by :meth:`PluginManager.process_parsers` for
        each plugin instance.

        Returns
        -------
        None
        """
        pass

    def help(self):
        """Return help for this plugin.

        The help string will be used in the command line option.
        """
        return "(no help available)"

    def preprocess(self, molecules, args):
        """Preprocess the molecules.

        Executed before process.

        Parameters
        ----------
        molecules: list of :obj:`mollib.Molecule` objects
            The list of molecule objects to pre-process.
        args: :obj:`argparse.ArgumentParser`
            The ArgumentParser parsed arguments.
        """
        pass

    def process(self, molecules, args):
        """Process the given molecules.

        Subclasses implement this method to conduct operations on the molecule.
        The parent method does not need to be called.

        Parameters
        ----------
        molecules: list of :obj:`mollib.Molecule` objects
            The list of molecule objects to pre-process.
        args: :obj:`argparse.ArgumentParser`
            The ArgumentParser parsed arguments.
        """
        pass

    def postprocess(self, molecules, args):
        """Postprocess the molecule.

        Executed after process.

        Parameters
        ----------
        molecules: list of :obj:`mollib.Molecule` objects
            The list of molecule objects to pre-process.
        args: :obj:`argparse.ArgumentParser`
            The ArgumentParser parsed arguments.
        """
        pass

