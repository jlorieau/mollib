"""
Analysis functions for observed and predicted RDCs/RACSs
"""
from math import sqrt
from collections import OrderedDict
from itertools import groupby

from numpy import std
import scipy.stats

from mollib.utils.ordered_set import OrderedSet
from mollib.utils.numbers import round_sig
from .utils import sort_key
from . import settings


#: add stats on Euler angles (and minimal degenerate)
def calc_summary(magnetic_interactions, Saupe_components, data, predicted):
    """Calculate the statistics between predicted and calculated RDCs and RACSs.

    Parameters
    ----------
    magnetic_interactions: list of dicts
        - A list of dicts, one for each molecule to be fit.
          See :class:`mollib.pa.process_molecule.Process`
    Saupe_components: dict
        See the output of :func:`mollib.pa.svd.calc_pa_SVD`
    data: dict
        - **key**: interaction labels (str)
        - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data values.
    predicted: dict
        - **key**: interaction labels (str)
        - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data values.

    Returns
    -------
    summary: :obj:`collections.OrderedDict`
        - 'Q': (float) the Q-factor of the fit
        - 'R': (float) the R-factor of the fit
        - 'RMS': (Hz/ppb) the root-mean square of the fit
    """
    # Prepare variables to collect statistics
    summary = OrderedDict()
    RSS = 0.         # Residual Sum Squared
    RSS_scaled = 0.  # Residual Sum Squared (scaled by DCC or RCSA)
    count = 0        # Count of the number of data points.

    # Calculate the overal Aa and Ar from the sum of each structural component
    Aa = Saupe_components['Aa']
    Ar = Saupe_components['Ar']
    sum_Aa = sum(Aa)
    sum_Ar = sum(Ar)
    sum_Rh = sum_Ar/sum_Aa

    # Loop over the data (observed values), and calculate the RSS with the
    # calculated values
    for key in data:
        if key not in predicted:
            continue
        obs = data[key].value
        calc = predicted[key].value
        value = next(i[key] for i in magnetic_interactions if key in i)
        if value is None:
            continue

        scale, _ = value

        residual = (obs - calc)**2
        RSS += residual
        RSS_scaled += residual / scale**2

        count += 1

    # Calculate the stats: Q-factor, R-factor, RSS.
    # Round these numbers to remove insignificant digits
    Q = 100. * sqrt(RSS_scaled / (float(count) * (sum_Aa)**2 *
                                  (4. + 3. * sum_Rh**2) / 5.))
    R = Q / sqrt(2.)
    RMS = sqrt(RSS / count)

    summary['Overall'] = OrderedDict()
    summary['Overall']['Q (%)'] = round(Q, 1)
    summary['Overall']['R (%)'] = round(R, 1)
    summary['Overall']['RSS'] = round(RSS, 1)
    summary['Overall']['RMS'] = round(RMS, 2)
    summary['Overall']['count']= count

    # Add statistics on each type interaction
    # Get the different interaction statistics
    sorted_keys = sorted(data.keys(), key=sort_key)
    interactions = OrderedSet()
    interactions.add('N-H')  # Add 'N-H' couplings default
    interactions |= [k[0] for k in map(sort_key, sorted_keys)]

    # Add basic stats for each type of interaction in the data
    for interaction in interactions:
        if interaction in settings.default_predicted_rdcs:
            summary[interaction] = OrderedDict()
            scale = settings.default_predicted_rdcs[interaction]
            summary[interaction]['Da (Hz)'] = round(scale * 2 * sum_Aa, 1)
            summary[interaction]['Rh'] = round(sum_Rh, 3)
        elif interaction in settings.default_predicted_racs:
            summary[interaction] = OrderedDict()
            scale = settings.default_predicted_rdcs[interaction]
            summary[interaction]['Da (ppm)'] = round(scale * sum_Aa, 1)
            summary[interaction]['Rh'] = round(sum_Rh, 3)
        else:
            continue

    # Add statistics on the Saupe matrix and round to 4 sig figs
    summary['Saupe'] = OrderedDict()
    summary['Saupe']['Aa'] = sum(Saupe_components['Aa'])
    summary['Saupe']['Ar'] = sum(Saupe_components['Ar'])
    summary['Saupe']['Szz'] = Saupe_components['Szz']
    summary['Saupe']['Syy'] = Saupe_components['Syy']
    summary['Saupe']['Sxx'] = Saupe_components['Sxx']

    for k,v in summary['Saupe'].items():
        summary['Saupe'][k] = round_sig(v, 4)

    # Add statistics on the Saupe Orientation and round to 1 decimal
    summary['Angles ZYZ (deg)'] = OrderedDict()
    summary['Angles ZYZ (deg)']['alpha'] = Saupe_components['alpha_z']
    summary['Angles ZYZ (deg)']['beta'] = Saupe_components['beta_y']
    summary['Angles ZYZ (deg)']['gamma'] = Saupe_components['gamma_z']

    for k,v in summary['Angles ZYZ (deg)'].items():
        summary['Angles ZYZ (deg)'][k] = round(v, 1)

    return summary


def G_critical(N, alpha=0.05):
    """Calculate the G-critical value (two-sided) for the given data count and
    alpha value.

    Parameters
    ----------
    N: int
        The number of data points.
    alpha: float
        The alpha critical value.

    Returns
    -------
    G: float
        The two-sided Grubbs test critical value.

    Examples
    --------
    >>> G = G_critical(10, 0.05)
    >>> round(G, 2)
    2.29
    >>> G = G_critical(17, 0.05)
    >>> round(G, 2)
    2.62
    >>> G = G_critical(40, 0.01)
    >>> round(G, 2)
    3.38
    """
    significance_level = alpha / (2. * N)
    t = scipy.stats.t.isf(significance_level, N - 2)
    return ((N - 1) / sqrt(N)) * (sqrt(t ** 2 / (N - 2 + t ** 2)))


def find_outliers(data, predicted):
    """Find data points that deviate from predicted values more than the others.

    Outliers are identified from a two-sided Grubbs test, using the alpha values
    (alpha_warning, alpha_bad) specified in the settings file. Outliers are
    only identified within each group of RDC or RACS.


    Returns
    -------
    (warning, bad): tuple of lists
        The lists of interaction labels for data points that are considered
        either warning outliers or bad outliers.
    """
    # Group the data and predicted by RDC or RACS type.
    # The sort_key will convert a label like '14N-H' into ('N-H', 14). To group
    # the data into rdc and racs types, we will use the first item of this
    # tuple to group the values.
    keys_sorted = sorted(data, key=lambda x: sort_key(x)[0])
    data_groups = {k:list(g) for k, g in
                   groupby(keys_sorted, key=lambda x: sort_key(x)[0])}

    # Keep track of the standard deviation for each group
    stdev_dict = {}
    count_dict = {}

    # Calculate the standard deviation for each group
    for group_name, labels in data_groups.items():
        values = []
        count = 0
        for label in sorted(labels, key=sort_key):
            if label in data and label in predicted:
                deviation = data[label].value - predicted[label].value
                values.append(deviation)
                count += 1
            else:
                continue

        # If no points were found for this group, then the stdev and count
        # cannot be calculated.
        if count <= 1:
            stdev_dict[group_name] = None
            count_dict[group_name] = None
        # Otherwise process the stdev and count
        else:
            stdev = std(values)
            stdev_dict[group_name] = stdev
            count_dict[group_name] = count

    # Go back and find all values that are outside of the 'warning' and bad
    # confidence interval
    warning = []
    bad = []
    for group_name, labels in data_groups.items():
        stdev = stdev_dict[group_name]
        count = count_dict[group_name]

        if stdev is None or count is None:
            continue

        # Calculate the Grubbs test critical cutoffs
        warning_cutoff = G_critical(count, settings.alpha_warning)
        bad_cutoff = G_critical(count, settings.alpha_bad)

        for label in sorted(labels, key=sort_key):
            if label not in data or label not in predicted:
                continue

            deviation = data[label].value - predicted[label].value

            # See if the deviation falls outside the 'warning' or 'bad'
            # G-critical values
            if abs(deviation/stdev) > bad_cutoff:
                bad.append(label)
            elif abs(deviation/stdev) > warning_cutoff:
                warning.append(label)
            else:
                continue

    return warning, bad
