"""
The base Plugin class.
"""

# Common arguments for all parsers
# in_args = ('-i',)
# in_kwargs = {'action': 'append',
#              'nargs': '+',
#              'required': True,
#              'type': str,
#              'metavar': 'identifier/filename',
#              'help': ('(required) The filename(s) or PDB identifier(s) of the'
#                       ' structure(s)')
#             }
#
# config_args = ('-c', '--config')
# config_kwargs = {'nargs': 1,
#                  'required': False,
#                  'type': str,
#                  'metavar': 'filename',
#                  'help': "The configuration file's name."
#                  }
#
# hydrogenate_args = ('--hydrogenate',)
# hydrogenate_kwargs = {'action': 'store_true',
#                       'help': 'Strip the hydrogens and re-add hydrogens'}
#

class Plugin(object):
    """Abstract base class for mollib plugins.

    Attributes
    ----------
    name: str
        The name of the plugin.
    enabled: bool
        If true, the plugin is enabled.
    command: str
        The shell command for this plugin.
    order: int
        The order to place the plugin in the command line options. (Higher is
        lower priority)
    parents: dict
        A dict of the parent processors. This is used if plugins want to add
        options to other commands. The key is the command (str) and the value
        is the parent parser. For example, the 'process' parent parser is
        contains the arguments for the 'process' command.


        .. note:: It is recommended that this dict be accessed with a
                  setdefault method and a default value of:
                  subparsers.add_parser('', add_help=False)


    .. note:: Preprocessor and postprocessor plugins should always return True
              for the selected method. These should also add to the parent
              parser.

    .. note:: Preprocessor and postprocessor plugins should have an order
              lower than 0 and should access the 'process' parent parser.
              Analysis methods should have an order between
              200-1000.
    """

    name = None
    enabled = True
    command = None
    order = 200
    parents = {}
    argument_title = 'arguments'

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
        if self.name is None:
            self.name = self.__class__.__name__.lower()
        if self.command is None:
            self.command = self.__class__.__name__.lower()

    def options(self, subparsers):
        """Register the command line options.

        The base method creates a parser from the passed subparsers, with name
        from self.command.
        The subclasses method should add arguments to this new parser.

        Parameters
        ----------
        parser: :obj:`argparse.ArgumentParser`
            A root argparse instance.

        Returns
        -------
        :obj:`argparse.ArgumentParser`
            The parser used by this plugin.
        """
        p = subparsers.add_parser(self.command, help=self.help(),
                                  parents=self.parents.values())
        p._optionals.title = self.argument_title

        return p

    def help(self):
        """Return help for this plugin.

        The help string will be used in the command line option.
        """
        return "(no help available)"

    def selected(self, args):
        """Test whether this plugin was selected.

        The application will run this plugin only if this method is True.

        Parameters
        ----------
        args: :obj:`argparse.ArgumentParser`
            The ArgumentParser parsed arguments.

        Returns
        -------
        bool
            True if this plugin was selected, otherwise False.
        """
        if (self.enabled and
            hasattr(args, 'command') and
            args.command == self.command):
            return True
        else:
            return False

    def preprocess(self, molecule, args):
        """Preprocess the molecule.

        Executed before process.

        Parameters
        ----------
        molecule: :obj:`mollib.Molecule`
            The molecule object to pre-process.
        args: :obj:`argparse.ArgumentParser`
            The ArgumentParser parsed arguments.
        """
        pass

    def process(self, molecule, args):
        """Process the given molecule.

        Subclasses implement this method to conduct operations on the molecule.
        The parent method does not need to be called.

        Parameters
        ----------
        molecule: :obj:`mollib.Molecule`
            The molecule object to process.
        args: :obj:`argparse.ArgumentParser`
            The ArgumentParser parsed arguments.
        """
        pass

    def postprocess(self, molecule, args):
        """Postprocess the molecule.

        Executed after process.

        Parameters
        ----------
        molecule: :obj:`mollib.Molecule`
            The molecule object to post-process.
        args: :obj:`argparse.ArgumentParser`
            The ArgumentParser parsed arguments.
        """
        pass

    @classmethod
    def plugin_instances(cls):
        if hasattr(cls, '_instances'):
            return sorted(cls._instances, key=lambda i: i.order)
        return []
