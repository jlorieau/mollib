class PrimitiveMetaClass(type):
    """A Primitive Metaclass for properly assigning and collecting attributes,
    such as 'optional'"""
    def __new__(mcs, classname, bases, class_dict):

        # This code adds 'optional' tuples from parent classes to the
        # 'optional' attribute of this class.
        parent_optional = tuple(*[getattr(base, 'optional')
                                  for base in bases
                                  if hasattr(base, 'optional')])
        class_dict['optional'] = class_dict['optional'] + parent_optional \
            if 'optional' in class_dict else parent_optional

        return type.__new__(mcs, classname, bases, class_dict)


class Primitive(object):
    """The base for all objects."""

    __metaclass__ = PrimitiveMetaClass

    __slots__ = ()
    optional = ()

    def __init__(self, **kwargs):

        # Check that the required arguments have been specified
        req_kwargs = [kw for kw in self.__slots__ if kw not in self.optional]
        assert all(kw in kwargs for kw in req_kwargs), \
            "All of the following parameters are needed: {}"\
            .format(req_kwargs)

        # Assign the values
        [setattr(self, kw, value) for kw, value in kwargs.items()
         if kw in self.__slots__]

        super(Primitive, self).__init__()
