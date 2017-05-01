# -*- coding: utf-8 -*-


class Datum(object):
    """A data value.

    Attributes
    ----------
    value: float
        The value of the data point.
    error: error
        The value of the error for the data point.
    """
    value = None
    error = None

    def __init__(self, **kwargs):
        for k,v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)

    def __repr__(self):
        if not hasattr(self, '_repr'):
            self._repr = "{}({:.1f}".format(self.__class__.__name__.lower(),
                                            self.value)
            if self.error is None:
                self._repr += ")"
            else:
                self._repr += "±{:.1f})".format(self.error)

        return self._repr
