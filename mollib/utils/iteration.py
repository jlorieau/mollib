try:
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest


def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx."
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)
