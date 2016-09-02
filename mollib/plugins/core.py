"""
The base Plugin class.
"""

# Common arguments for all parsers
in_args = ('-i',)
in_kwargs = {'action': 'append',
             'nargs': '+',
             'required': True,
             'type': str,
             'metavar': 'identifier/filename',
             'help': ('(required) The filename(s) or PDB identifier(s) of the'
                      ' structure(s)')
            }

argument_title = 'arguments'


class Plugin(object):
    """Base class for mollib plugins.

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
    """

    name = None
    enabled = True
    command = None
    order = 100

    def __new__(cls, *args, **kwargs):
        "Keep track of class instances"
        instance = object.__new__(cls, *args, **kwargs)
        if "_instances" not in cls.__dict__:
            cls._instances = []
        cls._instances.append(instance)
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
        parser: :obj:`argparse`
            A root argparse instance.

        Returns
        -------
        :obj:`argparse`
            The parser used by this plugin.
        """
        p = subparsers.add_parser(self.command, help=self.help())
        p._optionals.title = argument_title
        p.add_argument(*in_args, **in_kwargs)
        return p

    def help(self):
        """Return help for this plugin.

        The help string will be used in the command line option.
        """
        return "(no help available)"

    def process_molecule(self, molecule):
        """Process the given molecule.

        Subclasses implement this method to conduct operations on the molecule.
        The parent method does not need to be called.

        Parameters
        ----------
        molecule: :obj:`molecule`
            The molecule object to process.
        """
        raise NotImplementedError

    @classmethod
    def plugin_instances(cls):
        if hasattr(cls, '_instances'):
            return sorted(cls._instances, key=lambda i: i.order)
        return []
