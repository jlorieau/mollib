# -*- coding: utf-8 -*-
"""
Functions to conduct the SVD on the RDCs and RACS.
"""
import logging

import numpy as np
from scipy import linalg

from .analysis import calc_statistics
from .utils import get_data_type, sort_key
from . import settings


#: A set of labels for which the analysis has not been implemented. We keep
#: a set here so that the error message isn't printed multiple times by
#: calc_pa_SVD
not_implemented_errors = set()


def get_error(label, data):
    """Return the error for the given interaction.

    Parameters
    ----------
    label: str
        interaction label
    data: dict
        The experimental/observed RDC and RACS data.
        - **key**: interaction labels (str)
        - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data values.

    Returns
    -------
    error: float
        The interaction's error.
    """
    assert label in data

    # Use the data point's error, if it's specified. (i.e. it's not None or
    # equal to zero.)
    if data[label].error is not None and data[label].error != 0.0:
        return data[label].error

    # Otherwise calculate a default value
    interaction_type = sort_key(label)[0]
    if interaction_type in settings.default_error:
        value = data[label].value
        rel_error = settings.default_error[interaction_type]
        return value * rel_error

    interaction_type_rev = '-'.join(interaction_type.split('-')[::-1])
    if interaction_type_rev in settings.default_error:
        value = data[label].value
        rel_error = settings.default_error[interaction_type_rev]
        return value * rel_error

    msg = "Error of type '{}' not found for '{}'"
    logging.error(msg.format(interaction_type, label))

    return 1.0


# TODO: incorporate errors
def calc_pa_SVD(magnetic_interactions, data):
    """Calculate the best-fit Saupe matrices for the given magnetic interaction
    arrays and RDC/RACS data.

    Parameters
    ----------
    magnetic_interactions: list of dicts
        A list of dicts, one for each molecule to be fit.
        See :class:`mollib.pa.process_molecule.Process`.

    data: dict
        The experimental/observed RDC and RACS data.
        - **key**: interaction labels (str)
        - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data values.

    Returns
    -------
    (data_pred, Saupe_components, stats): tuple
        - data_pred: dict
              The SVD predicted RDC and RACS data.
              - **key**: interaction labels (str)
              - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data
                values.
        - Saupe_components: dict
              - 'S_xyz': (list of arrays) The 3x1 Saupe matrix in x/y/z repr.
              - 'Aa': (list) The degree of alignment.
              - 'Ar': (list) The alignment rhombicity.
              - 'Rh': (list) The rhombicity
        - stats: dict
              - See :func:`calc_statistics`
    """
    assert(isinstance(magnetic_interactions, list))
    global not_implemented_errors

    # Make an ordered list of the keys
    ordered_keys = sorted(data.keys())

    # Construct the A and D matrices
    A = []
    D = []
    for key in ordered_keys:
        for interaction_dict in magnetic_interactions:
            # If the key isn't in the interaction_dict, then it's not known
            # how to process this interaction
            if key not in interaction_dict:
                if key not in not_implemented_errors:
                    msg = "Processing of data point '{}' is not implemented."
                    logging.error(msg.format(key))
                    not_implemented_errors.add(key)
                continue

            scale, arr = interaction_dict[key]
            A.extend([arr,])

            expt_value = data[key].value
            expt_error = get_error(key, data)
            print(key, data[key], data[key].error == 0.0)
            D.append(expt_value)


    # Create an array from the A and D matrices
    A = np.array(A)
    D = np.array(D)

    # conduct the SVD on the A matrix
    U, w, V = linalg.svd(A, full_matrices=False)

    # Remove inf items from the inverted 'w' vector
    w_inv = 1./w
    for i in range(w_inv.size):
        w_inv[i] = 0. if np.isinf(w_inv[i]) else w_inv[i]
    w_inv = linalg.diagsvd(w_inv, w.shape[0], w.shape[0])

    # Calculate the S matrix
    A_inv = np.dot(V.transpose(), np.dot(w_inv, U.transpose()))
    S = np.dot(A_inv, D)

    # Calculate the A matrix for all the interactions in the interaction_dict.
    keys = {i for d in magnetic_interactions for i in d.keys()}
    ordered_keys = sorted(keys)
    A = []
    for key in ordered_keys:
        for interaction_dict in magnetic_interactions:
            if key not in interaction_dict:
                continue
            scale, arr = interaction_dict[key]
            A.extend([arr,])

    D_pred = np.dot(A, S)

    # Calculate the predicted data. It must be sorted the same way as the
    # interaction_dict
    data_pred = {}
    for key, D in zip(ordered_keys, D_pred):
        # Determine whether the predicted data is an RDC or RACS
        data_type = get_data_type(key)
        data_pred[key] = data_type(value=D, error=0.0)

    # Break up the S matrix into individual Saupe matrices, Das and Rh
    Saupe_components = {}
    for i in ('S_xyz', 'Aa', 'Ar', 'Rh'):
        Saupe_components[i] = []

    for x in range(0, len(S), 5):
        s = S[x:x + 5]
        s_xyz = np.array([[-0.5 * (s[0] - s[1]), s[2], s[3], ],
                         [s[2], -0.5 * (s[0] + s[1]), s[4]],
                         [s[3], s[4], s[0]]])
        s_xyz = linalg.eigvals(s_xyz).real
        Saupe_components['S_xyz'].append(s_xyz)

        xx, yy, zz = [i for i in sorted(abs(s_xyz))]
        aa = max(s_xyz) / 2. if max(s_xyz) == zz else min(s_xyz) / 2.
        ar = (yy - xx) / 3.
        Saupe_components['Aa'].append(aa)
        Saupe_components['Ar'].append(ar)
        Saupe_components['Rh'].append(ar / abs(aa))

    # Calculate the statistics
    stats = calc_statistics(magnetic_interactions, Saupe_components, data,
                            data_pred)

    return (data_pred, Saupe_components, stats)


