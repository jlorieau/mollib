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
from mollib.utils.interactions import sort_func, interaction_type
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
        - **'Overall'**: Overall Statistics
            - 'Q (%)': (float) The fit Q-factor in percentage
            - 'RMS': (Hz/ppb) The root-mean square of the fit
            - 'count': (int) The number of interactions fit
        - **'Alignment'**: Details on the alignment tensor
            - 'Aa': (float) The alignment tensor anisotropy
            - 'Ar': (float) The alignment tensor rhobicity
        - **'Saupe'**: Details on the Saupe matrix
            - 'Szz': (float) The zz-component of the Saupe matrix
            - 'Sxx': (float) The xx-component of the Saupe matrix
            - 'Syy': (float) The yy-component of the Saupe matrix
        - **'Angles'**: Alignment tensor orientation in Rose convention
            - "Z (deg)": (degrees) The alignment alpha angle
            - "Y' (deg)": (degrees) The alignment beta angle
            - "Z'' (deg)": (degrees) The alignment gamma angle
    """
    # Calculate the overal Aa and Ar from the sum of each structural component
    Aa = Saupe_components['Aa']
    Ar = Saupe_components['Ar']
    sum_Aa = sum(Aa)
    sum_Ar = sum(Ar)
    sum_Rh = sum_Ar/sum_Aa

    # Prepare variables to collect statistics
    summary = OrderedDict()
    RSS = {}  # Residual Sum Squared
    RSS_scaled = {}  # Residual Sum Squared (scaled by DCC or RCSA)
    count = {}  # Count of the number of data points.

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

        # Get the value from the data
        scale, _ = value

        # Identify the interaction type (str) for this data value
        label = interaction_type(key)

        # Calculate the overall and interaction_type specific statistics
        residual = (obs - calc)**2

        # Calculate the overall RSS and interaction-specific RSS
        RSS['Overall'] = RSS.setdefault('Overall', 0.0) + residual
        RSS[label] = RSS.setdefault(label, 0.0) + residual

        # Calculate the overall RSS_scaled and interaction-specific RSS
        RSS_scaled['Overall'] = (RSS_scaled.setdefault('Overall', 0.0) +
                                 residual / scale**2)
        RSS_scaled[label] = (RSS_scaled.setdefault(label, 0.0) +
                                 residual / scale**2)

        # Add it to the overall and interaction-specific counts
        count['Overall'] = count.setdefault('Overall', 0) + 1
        count[label] = count.setdefault(label, 0) + 1

    # Calculate the Overall and interaction-specific Q-factor and RMS.
    Q = {}
    RMS = {}
    for key in RSS:
        Q[key] = 100. * sqrt(RSS_scaled[key] / (float(count[key]) *
                                                (sum_Aa)**2 *
                                                (4. + 3. * sum_Rh**2) / 5.))
        RMS[key] = sqrt(RSS[key] / float(count[key] - 1))

    # Round these numbers to remove insignificant digits and add it to the
    # summary (result) dict
    summary['Overall'] = OrderedDict()
    summary['Overall']['Q (%)'] = round(Q['Overall'], 1)
    summary['Overall']['RMS'] = round(RMS['Overall'], 2)
    summary['Overall']['count']= count['Overall']

    # Add statistics on each type interaction
    # Get the different interaction statistics
    sorted_keys = sorted(data.keys(), key=sort_func)
    interactions = OrderedSet()
    interactions.add('N-H')  # Add 'N-H' couplings default
    interactions |= [k for k in map(interaction_type, sorted_keys)]

    # Add basic stats for each type of interaction in the data
    for interaction in interactions:
        summary[interaction] = OrderedDict()

        if interaction in Q:
            summary[interaction]['Q (%)'] = round(Q[interaction], 1)
        if interaction in RMS:
            summary[interaction]['RMS'] = round(RMS[interaction], 2)
        if interaction in count:
            summary[interaction]['count'] = count[interaction]

        # Now add the Da/Rh for each interaction type
        if interaction in settings.default_predicted_rdcs:
            scale = settings.default_predicted_rdcs[interaction]
            summary[interaction]['Da (Hz)'] = round(scale * 2 * sum_Aa, 1)
            summary[interaction]['Rh'] = round(sum_Rh, 3)
        elif interaction in settings.default_predicted_racs:
            scale = settings.default_predicted_racs[interaction]['delta']
            summary[interaction]['Da (ppb)'] = round(scale * 1000 * sum_Aa, 1)
            summary[interaction]['Rh'] = round(sum_Rh, 3)
        else:
            continue

    # Add statistics on the Saupe matrix and round to 4 sig figs
    summary['Alignment'] = OrderedDict()
    summary['Alignment']['Aa'] = sum(Saupe_components['Aa'])
    summary['Alignment']['Ar'] = sum(Saupe_components['Ar'])

    for k,v in summary['Alignment'].items():
        summary['Alignment'][k] = round_sig(v, 4)

    summary['Saupe'] = OrderedDict()
    summary['Saupe']['Szz'] = Saupe_components['Szz']
    summary['Saupe']['Syy'] = Saupe_components['Syy']
    summary['Saupe']['Sxx'] = Saupe_components['Sxx']

    for k,v in summary['Saupe'].items():
        summary['Saupe'][k] = round_sig(v, 4)

    # Add statistics on the Saupe Orientation and round to 1 decimal
    summary['Angles'] = OrderedDict()
    summary['Angles']["Z (deg)"] = Saupe_components['alpha_z']
    summary['Angles']["Y' (deg)"] = Saupe_components['beta_y']
    summary['Angles']["Z'' (deg)"] = Saupe_components['gamma_z']

    for k,v in summary['Angles'].items():
        summary['Angles'][k] = round(v, 1)

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
    # The sort_func will convert a label like '14N-H' into ('N-H', 14). To group
    # the data into rdc and racs types, we will use the first item of this
    # tuple to group the values.
    keys_sorted = sorted(data, key=lambda x: sort_func(x)[0])
    data_groups = {k:list(g) for k, g in
                   groupby(keys_sorted, key=lambda x: sort_func(x)[0])}

    # Keep track of the standard deviation for each group
    stdev_dict = {}
    count_dict = {}

    # Calculate the standard deviation for each group
    for group_name, labels in data_groups.items():
        values = []
        count = 0
        for label in sorted(labels, key=sort_func):
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

        for label in sorted(labels, key=sort_func):
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
