
def identify_hydridization(atom):
    """The hybridization string for the given atom.

    Parameters
    ----------
    atom : :obj:`atom`
        The atom to identify the hydridization for.

    Returns
    -------
    str
        The string for the hybridization of this atom. ex: 'sp3', 'sp2'
    """
    if atom.element in ('C', 'N', 'O'):
        number_bonded = len(atom.topology)

        return 'sp' + (str(number_bonded - 1) if number_bonded > 1 else '')
    if atom.element == 'H':
        return 's'
