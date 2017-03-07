from .exceptions import ParameterError


msg = ("Either 3 tensor component arguments, the delta/eta kwargs or "
        "span/skew kwargs should be given.")


def get_Haeberlen(*args, **kwargs):
    """Calculate the tensor components and values using the Haeberlen
    convention.

        | dzz - iso | >= | dxx - iso | >= | dyy - iso |

        delta = dzz - iso

        eta = (dxx - dyy) / delta

        iso = (dzz + dxx + dyy) / 3.

    Parameters
    ----------
    This function can be used one of three ways:
    args:
        1. The 3 components of the tensor (xx/yy/zz or 11/22/33), in any order.
    kwargs:
        2. The delta and eta kwargs. The iso kwarg is optional .
        3. The span and skew kwargs. The iso kwarg is optional.

    Returns
    -------
    (delta, eta, iso, dzz, dxx, dyy): tuple of floats

    Examples
    --------
    >>> get_Haeberlen(74.7, 11.8, -86.5)
    (-86.5, 0.727, 0.0, -86.5, 74.7, 11.8)
    >>> get_Haeberlen(delta=-86.5, eta=0.727)
    (-86.5, 0.727, 0.0, -86.5, 74.69, 11.81)
    >>> get_Haeberlen(span=161.2, skew=0.22)
    (-86.51, 0.727, -0.0, -86.51, 74.69, 11.82)
    >>> get_Haeberlen(108.5, -45.7, -62.8)
    (108.5, 0.158, 0.0, 108.5, -62.8, -45.7)
    >>> get_Haeberlen(delta=108.5, eta=0.158)
    (108.5, 0.158, 0.0, 108.5, -62.82, -45.68)
    >>> get_Haeberlen(span=171.32, skew=-0.8)
    (108.5, 0.158, -0.0, 108.5, -62.82, -45.69)
    """
    if len(args) == 3:
        pass
    elif 'delta' in kwargs and 'eta' in kwargs:
        args = get_Herzfeld(**kwargs)[3:]
    elif 'span' in kwargs and 'skew' in kwargs:
        args = get_Herzfeld(**kwargs)[3:]
    else:
        global msg
        raise ParameterError(msg)

    iso = round(sum(args) / len(args), 2)
    dyy, dxx, dzz = sorted(args, key=lambda x: abs(x - iso))
    delta = dzz - iso
    eta = round(abs((dxx - dyy) / delta), 3)

    return delta, eta, iso, dzz, dxx, dyy


def get_Herzfeld(*args, **kwargs):
    """Calculate the tensor components and values using the Herzfeld
    convention.

        d11  >=  d22  >=  d33

        span = d11 - d33

        skew = 3. * (d22 - iso) / span

        iso = (d11 + d22 + d33) / 3.

    Parameters
    ----------
    This function can be used one of three ways:
    args:
        1. The 3 components of the tensor (xx/yy/zz or 11/22/33), in any order.
    kwargs:
        2. The delta and eta kwargs. The iso kwarg is optional.
        3. The span and skew kwargs. The iso kwarg is optional.

    Returns
    -------
    (span, skew, iso, d11, d22, d33): tuple of floats

    Examples
    --------
    >>> get_Herzfeld(74.7, 11.8, -86.5)
    (161.2, 0.22, 0.0, 74.7, 11.8, -86.5)
    >>> get_Herzfeld(delta=-86.5, eta=0.727)
    (161.19, 0.22, 0.0, 74.69, 11.81, -86.5)
    >>> get_Herzfeld(span=161.2, skew=0.22)
    (161.2, 0.22, 0.0, 74.69, 11.82, -86.51)
    >>> get_Herzfeld(108.5, -45.7, -62.8)
    (171.3, -0.8, 0.0, 108.5, -45.7, -62.8)
    >>> get_Herzfeld(delta=108.5, eta=0.158)
    (171.32, -0.8, 0.0, 108.5, -45.68, -62.82)
    >>> get_Herzfeld(span=171.32, skew=-0.8)
    (171.32, -0.8, 0.0, 108.5, -45.69, -62.82)
    """
    if len(args) == 3:
        iso = sum(args) / len(args)
        d33, d22, d11 = sorted(args, key=lambda x: x - iso)

    elif 'delta' in kwargs and 'eta' in kwargs:
        delta = kwargs['delta']
        eta = kwargs['eta']
        iso = kwargs.get('iso', 0.0)
        if delta > 0:
            d11 = iso + delta
            d22 = iso - delta * (1. - eta) / 2.
            d33 = iso - delta * (1. + eta) / 2.
        else:
            d33 = iso + delta
            d22 = iso - delta * (1. - eta) / 2.
            d11 = iso - delta * (1. + eta) / 2.

        d11 = round(d11, 2)
        d22 = round(d22, 2)
        d33 = round(d33, 2)

    elif 'span' in kwargs and 'skew' in kwargs:
        span = kwargs['span']
        skew = kwargs['skew'] * -1.
        iso = kwargs.get('iso', 0.0)

        d22 = iso - span * skew / 3.
        d33 = (3. * iso - d22 - span) / 2.
        d11 = 3. * iso - d22 - d33

        d11 = round(d11, 2)
        d22 = round(d22, 2)
        d33 = round(d33, 2)

    else:
        global msg
        raise ParameterError(msg)

    span = d11 - d33
    skew = round(3. * (d22 - iso) / span, 3)

    return span, skew, iso, d11, d22, d33
