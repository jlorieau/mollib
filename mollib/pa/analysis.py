"""
Analysis functions for observed and predicted RDCs/RACSs
"""
from math import sqrt
from collections import OrderedDict

# class Summary(object):
#     pass
#
# class Results(Summary):
#     pass


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


    # Calculate the stats: Q-factor, R-factor, RSS.
    # Round these numbers to remove insignificant digits
    stats['Q'] = sqrt(RSS_scaled /
                     (float(count) * (sum_Aa)**2 * (4. + 3. * sum_Rh**2) /5.))
    stats['R'] = stats['Q'] / sqrt(2.)
    stats['RSS'] = RSS

    for k,v in stats.items():
        stats[k] = round(v, 3)

    stats['count'] = count

    return stats