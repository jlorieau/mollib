from .utils import sort_key
from .analysis import find_outliers
from mollib.utils import MDTable, FormattedStr


def center(number, rjust=4):
    """Center a number string about the decimal."""
    number = str(number)
    split = number.split('.')
    return (''.join((split[0].rjust(4), '.', split[1]))
            if len(split) > 1 else split[0])


def report_tables(data, predicted=None):
    """Produce the partial alignment report for the observed and predicted
    RDC and RACS values.

    Parameters
    ----------
    data: dict
        The experimental/observed RDC and RACS data.
        - **key**: interaction labels (str)
        - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data values.
    predicted: dict (optional)
        The SVD predicted RDC and RACS data.
        - **key**: interaction labels (str)
        - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` predicted
          values.

    Returns
    -------
    tables: dict
        A dict with the tables:
            - **keys**
                - 'fit': the fit table
                - 'xx_pred': the predicted data table. ex: 'N-H_pred',
                  'CA-HA_pred'
            - **values**
                - The `mollib.utils.markdown.MDTable` objects.
    """
    if predicted is None:
        predicted = dict()

    # Prepare the fit data table
    tables = {}
    tables['fit'] = MDTable('Interaction', 'Value', 'Error',
                            'Predicted', 'Deviation')

    # Make a (shallow) copy of the predicted data dict. As we print off
    # interactions from the data, we remove them from this copied dict so that
    # we do not print the RDC/RACS twice.
    predicted_copy = predicted.copy()

    # Find the warning and bad outliers
    warning, bad = find_outliers(data, predicted)

    # Iterate over the data and add the values to the table.
    for label in sorted(data, key=sort_key):
        # Get the fields
        interaction = label
        value = data[label].value
        error = data[label].error if data[label].error else '-'


        # Find the number of digits in the observed value so that the predicted
        # values (and deviations) can be rounded to the same precision
        split = str(value).split('.')
        if len(split) > 1:
            no_digits = len(str(value).split('.')[1])
        else:
            no_digits = 0

        # Get the predicted value and deviation between observed and predicted
        # (or put a '-' if there is None).
        if label in predicted_copy:
            pred = predicted_copy.pop(label).value
            pred = round(pred , no_digits)
            deviation = round(abs(value - pred), no_digits)
        else:
            pred = "-"
            deviation = "-"

        # Identify outlier points with either a warning (yellow) or as bad (red)
        # Also put an asterisk or exclamation mark next to the label's name
        if label in warning:
            fmt = 'yellow'
            interaction = label + '*'
        elif label in bad:
            fmt = 'red'
            interaction = label + '!'
        else:
            fmt = ''

        # Add the fields to a row in the table
        tables['fit'].add_row(FormattedStr(interaction, fmt),
                              FormattedStr(center(value), fmt),
                              FormattedStr(center(error), fmt),
                              FormattedStr(center(pred), fmt),
                              FormattedStr(center(deviation), fmt))

    # Prepare tables for predicted values
    predicted_interactions = set([sort_key(i)[2]
                                  for i in predicted_copy.keys()])

    # Populate the table and rowsrows
    for label in sorted(predicted_copy.keys(), key=sort_key):
        # Get the appropriate table
        interaction_type = sort_key(label)[0]
        table = tables.setdefault(interaction_type + '_pred',
                                  MDTable('Interaction', 'Predicted'))

        # Get the fields
        value = round(predicted_copy[label].value, 2)

        # Add the fields to a row in the table
        table.add_row(label, center(value))

    return tables