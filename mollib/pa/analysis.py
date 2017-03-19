"""
Analysis functions for observed and predicted RDCs/RACSs
"""
from math import sqrt
from collections import OrderedDict
from itertools import groupby

from numpy import std
from scipy import stats

from .utils import sort_key
from . import settings


# TODO: Add groupwise statistics, including average Das
def calc_statistics(magnetic_interactions, Saupe_components, data, predicted):
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
    stats: :obj:`collections.OrderedDict`
        - 'Q': (float) the Q-factor of the fit
        - 'R': (float) the R-factor of the fit
        - 'RMS': (Hz/ppb) the root-mean square of the fit
    """
    # Prepare variables to collect statistics
    stats = OrderedDict()
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

    # for rdc_type, static_value in settings.default_predicted_rdcs.items():
    #     if rdc_type not in settings.default_error:
    #         continue
    #     error = settings.default_error[rdc_type]
    #     print(rdc_type, static_value, static_value * 2 * sum_Aa, sum_Rh )

    # Add statistics on the Saupe matrix
    stats['Da H-N (Hz)'] = settings.default_predicted_rdcs['N-H'] * 2. * sum_Aa
    stats['Da H-N (Hz)'] = round(stats['Da H-N (Hz)'], 2)
    stats['Rh'] = sum_Rh
    stats['Rh'] = round(stats['Rh'], 3)

    # Calculate the stats: Q-factor, R-factor, RSS.
    # Round these numbers to remove insignificant digits
    stats['Q-factor (%)'] = 100. *sqrt(RSS_scaled /
                            (float(count) * (sum_Aa)**2 *
                             (4. + 3. * sum_Rh**2) /5.))
    stats['Q-factor (%)'] = round(stats['Q-factor (%)'], 1)

    stats['R-factor (%)'] = stats['Q-factor (%)'] / sqrt(2.)
    stats['R-factor (%)'] = round(stats['R-factor (%)'], 1)

    stats['RSS'] = RSS
    stats['RSS'] = round(stats['RSS'], 1)

    stats['RMS'] = sqrt(RSS / count )
    stats['RMS'] = round(stats['RMS'], 2)

    stats['count'] = count

    return stats


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
    t = stats.t.isf(significance_level, N - 2)
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
