from itertools import groupby

from .utils import sort_key
from .analysis import find_outliers
from mollib.utils import MDTable, FormattedStr


def match_digits(ref, target):
    """Match the number of digits in a target number to """


def report_tables(data, predicted):
    """Produce the partial alignment report for the observed and predicted
    RDC and RACS values.

    Parameters
    ----------
    data: dict
        The experimental/observed RDC and RACS data.
        - **key**: interaction labels (str)
        - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data values.
    predicted: dict
        The SVD predicted RDC and RACS data.
        - **key**: interaction labels (str)
        - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` predicted
          values.

    Returns
    -------

    """
    # Group the data and predicted by RDC or RACS type.
    # The sort_key will convert a label like '14N-H' into ('N-H', 14). To group
    # the data into rdc and racs types, we will use the first item of this
    # tuple to group the values.
    data_groups = [k for k,g in
                   groupby(data, key=lambda x: sort_key(x)[0])]

    # Likewise, do the same grouping for the predicted_data
    # predicted_groups = [k for k,g in
    #                     groupby(predicted, key=lambda x: sort_key(x)[0])]
    #
    # Prepare the fit data table
    table_fit = MDTable('Interaction', 'Value', 'Error',
                         'Predicted', 'Deviation')

    # Make a (shallow) copy of the predicted data dict. As we print off
    # interactions from the data, we remove them from this copied dict so that
    # we do not print the RDC/RACS twice.
    predicted_copy = predicted.copy()

    # Find the warning and bad outliers
    warning, bad = find_outliers(data, predicted)

    # Iterate over the data and add the values to the table.
    for label in sorted(data, key=sort_key):
        interaction = label
        value = data[label].value
        error = data[label].error if data[label].error else '-'

        if label in predicted_copy:
            pred = predicted.pop(label).value
            deviation = abs(value - pred)
        else:
            pred = "-"
            deviation = "-"

        if label in warning:
            fmt = 'yellow'
        elif label in bad:
            fmt = 'red'
        else:
            fmt = ''

        table_fit.add_row(FormattedStr(interaction, fmt),
                          FormattedStr(value, fmt),
                          FormattedStr(error, fmt),
                          FormattedStr(pred, fmt),
                          FormattedStr(deviation, fmt))

    print(table_fit.content())


