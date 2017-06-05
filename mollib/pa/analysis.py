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


# TODO: This function should be split into smaller functions, possibly a
# chain-of-command object pattern
def calc_summary(magnetic_interactions, Saupe_components, data, predicted):
    """Calculate the statistics between predicted and calculated RDCs and
    RACSs.

    Parameters
    ----------
    magnetic_interactions: list of dicts
        - A list of dicts, one for each molecule to be fit.
          See :class:`mollib.pa.process_molecule.Process`
    Saupe_components: dict
        See the output of :func:`mollib.pa.svd.calc_pa_SVD`
    data: dict
        - **key**: interaction labels, str
        - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data
          values.
    predicted: dict
        - **key**: interaction labels, str
        - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data
          values.

    Returns
    -------
    summary: :obj:`collections.OrderedDict`

        - 'Overall': Overall Statistics, :obj:`collections.OrderedDict`

          - 'Q (%)': The fit Q-factor in percentage, float
          - 'RMS': The root-mean square of the fit (Hz/ppb), float
          - 'count': The number of interactions fit, int

        - 'Alignment': Details on the alignment tensor,
          :obj:`collections.OrderedDict`

          - 'Aa': The alignment tensor anisotropy, float
          - 'Ar': The alignment tensor rhobicity, float

        - 'Saupe': Details on the Saupe matrix, :obj:`collections.OrderedDict`

          - 'Szz': The zz-component of the Saupe matrix, float
          - 'Sxx': The xx-component of the Saupe matrix, float
          - 'Syy': The yy-component of the Saupe matrix, float

        - 'Angles': Alignment tensor orientation in Rose convention,
          :obj:`collections.OrderedDict`

          - "Z (deg)": The alignment alpha angle (deg), float
          - "Y' (deg)": The alignment beta angle (deg), float
          - "Z'' (deg)": The alignment gamma angle (deg), float
    """
    # Calculate the overal Aa and Ar from the sum of each molecule/conformer
    Aa = [v for k, v in Saupe_components.items() if k.startswith('Aa')]
    Ar = [v for k, v in Saupe_components.items() if k.startswith('Ar')]

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
        Aarray = next(i[key] for i in magnetic_interactions if key in i)
        if Aarray is None:
            continue

        # Get the value from the data
        scale, _ = Aarray

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
                                                sum_Aa ** 2 *
                                                (4. + 3. * sum_Rh**2) / 5.))

        if count[key] > 1:
            RMS[key] = sqrt(RSS[key] / float(count[key] - 1))
        else:
            RMS[key] = '-'

    # Round these numbers to remove insignificant digits and add it to the
    # summary (result) dict
    summary['Overall'] = OrderedDict()
    summary['Overall']['Q (%)'] = round(Q['Overall'], 1)
    summary['Overall']['RMS'] = round(RMS['Overall'], 2)
    summary['Overall']['count'] = count['Overall']

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
            summary[interaction]['RMS'] = (round(RMS[interaction], 2)
                                           if isinstance(RMS[interaction],
                                                         float)
                                           else '-')
        if interaction in count:
            summary[interaction]['count'] = count[interaction]

    # Add statistics on the Saupe matrix and round to 4 sig figs.
    summary['Alignment'] = OrderedDict()
    summary['Alignment']['Aa'] = sum_Aa
    summary['Alignment']['Ar'] = sum_Ar

    for k,v in summary['Alignment'].items():
        summary['Alignment'][k] = round_sig(v, 4)

    # Put in statistics on the Saupe matrix and the alignment Da/Rh.
    # There is (potentially) one entry for each molecule/conformer identifier.
    no_molecules = len([i for i in Saupe_components.keys()
                        if i.startswith('Szz')])

    if no_molecules > 1:
        # In this case, there are multiple molecules conformers. An entry for
        # each must be created. Get the molecular identifiers for each,
        # extracting. These are, for example, ' (1)', ' (2)' and so on.
        ids = [i[3:] for i in Saupe_components.keys()
               if i.startswith('Szz')]
    else:
        # Otherwise there is just one molecule. The id for the molecule is
        # just an empty string
        ids = ['',]

    # Now add the Da/Rh for each interaction type
    for interaction in interactions:
        for id_ in sorted(ids):
            if interaction + id_ not in summary:
                summary[interaction + id_] = OrderedDict()

            Aa = Saupe_components['Aa' + id_]
            Rh = Saupe_components['Rh' + id_]

            if interaction in settings.default_predicted_rdcs:
                scale = settings.default_predicted_rdcs[interaction]
                Da = scale * 2 * Aa

                summary[interaction + id_]['Da' + id_ + ' (Hz)'] = round(Da, 1)
                summary[interaction + id_]['Rh'] = round(Rh, 3)

            elif interaction in settings.default_predicted_racs:
                scale = settings.default_predicted_racs[interaction]['delta']
                Da = scale * Aa * 1000.  # Convert from ppm to ppb

                summary[interaction + id_]['Da' + id_ + ' (ppb)'] = round(Da, 1)
                summary[interaction + id_]['Rh'] = round(Rh, 3)
            else:
                continue

    # Prepare statistics on the Saupe matrices
    for id_ in sorted(ids):
        summary['Saupe' + id_] = OrderedDict()
        summary['Saupe' + id_]['Szz'] = Saupe_components['Szz' + id_]
        summary['Saupe' + id_]['Syy'] = Saupe_components['Syy' + id_]
        summary['Saupe' + id_]['Sxx'] = Saupe_components['Sxx' + id_]

        # Round the numbers
        for k, v in summary['Saupe' + id_].items():
            summary['Saupe' + id_][k] = round_sig(v, 4)

    # Prepare statistics on the Saupe matrix angles
    for id_ in sorted(ids):
        summary['Angles' + id_] = OrderedDict()
        summary['Angles' + id_]["Z (deg)"] = Saupe_components['alpha_z'
                                                              + id_]
        summary['Angles' + id_]["Y' (deg)"] = Saupe_components['beta_y'
                                                               + id_]
        summary['Angles' + id_]["Z'' (deg)"] = Saupe_components['gamma_z'
                                                                + id_]

        for k, v in summary['Angles' + id_].items():
            summary['Angles' + id_][k] = round(v, 1)

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
    data_groups = OrderedDict((k, list(g)) for k, g in
                              groupby(keys_sorted,
                                      key=lambda x: sort_func(x)[0]))

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
