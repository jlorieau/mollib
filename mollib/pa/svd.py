# -*- coding: utf-8 -*-
"""
Functions to conduct the SVD on the RDCs and RACS.
"""
import logging

import numpy as np
from scipy import linalg

from mollib.utils.rotations import euler_zyz
from mollib.utils.interactions import sort_func
from .analysis import calc_summary
from .utils import get_data_type
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
    # Use the data point's error, if it's specified. (i.e. it's not None or
    # equal to zero.)
    if (label in data and data[label].error is not None and
        data[label].error != 0.0):
        return data[label].error

    # Otherwise calculate a default value
    interaction_type = sort_func(label)[0]
    if interaction_type in settings.default_error:
        return settings.default_error[interaction_type]

    interaction_type_rev = '-'.join(interaction_type.split('-')[::-1])
    if interaction_type_rev in settings.default_error:
        return settings.default_error[interaction_type_rev]

    msg = ("The default error for the interaction type '{}' was not "
           "found for '{}'")
    logging.info(msg.format(interaction_type, label))

    return 1.0


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
              - 'alpha_z': (float) The ZYZ alpha angle (degrees) of the Saupe
                 matrix
              - 'beta_y': (float) The ZYZ beta angle (degrees) of the Saupe
                 matrix
              - 'gamma_z': (float) The ZYZ gamma angle (degrees) of the Saupe
                 matrix
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
        # # If the key isn't in the interaction_dict, then it's not known
        # # how to process this interaction
        # if key not in interaction_dict:
        #     if key not in not_implemented_errors:
        #         msg = "Processing of data point '{}' is not implemented."
        #         logging.error(msg.format(key))
        #         not_implemented_errors.add(key)
        #     continue

        # Get the experimental value and error
        expt_value = data[key].value
        expt_error = get_error(key, data)

        D.append(expt_value / expt_error)

        # Construct the A-matrix
        A_line = []
        for interaction_dict in magnetic_interactions:
            # Check to see if the interaction has been processed.
            if key not in interaction_dict:
                if key not in not_implemented_errors:
                        msg = ("Processing of data point '{}' is not "
                               "implemented.")
                        logging.error(msg.format(key))
                        not_implemented_errors.add(key)
                continue

            scale, arr = interaction_dict[key]
            A_line.extend(arr * scale / expt_error)

        A.append(A_line)

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

    # Calculate the predicted dipolar couplings. It must be sorted the same way
    # as the D and A matrices
    D_pred = np.dot(A, S)
    data_pred = {}

    for key, D in zip(ordered_keys, D_pred):
        # Determine whether the predicted data is an RDC or RACS
        data_type = get_data_type(key)
        expt_error = get_error(key, data)
        data_pred[key] = data_type(value=D * expt_error, error=expt_error)

    # Break up the S matrix into individual Saupe matrices, Das and Rh. There
    # is one Saupe matrix for each molecule
    if len(magnetic_interactions) > 1:
        # This is the case for many molecules being fit
        molecule_number = range(1, len(magnetic_interactions) + 1)
    else:
        # This is the case for one molecule being fit.
        molecule_number = [None,]
    
    Saupe_components = {}
    if len(magnetic_interactions) > 1:

    else:
        for i in ('S_xyz', 'Aa', 'Ar', 'Rh'):
            Saupe_components[i] = []

    for x in range(0, len(S), 5):
        s = S[x:x + 5]

        # The order of the Saupe matrix components are:
        # Cyy, Czz, Cxy, Cxz, Cyz
        s_xyz = np.array([[-s[0]-s[1], s[2], s[3]],
                          [s[2],         s[0], s[4]],
                          [s[3],         s[4], s[1]]])

        # Conduct the eigen value decomposition
        (eig_values, eig_vectors) = linalg.eig(s_xyz)

        # Get the Saupe matrix angles
        alpha, beta, gamma = euler_zyz(eig_vectors)
        Saupe_components['alpha_z'] = alpha
        Saupe_components['beta_y'] = beta
        Saupe_components['gamma_z'] = gamma

        # Get the Saupe matrix values
        s_xyz = eig_values.real
        Saupe_components['S_xyz'].append(s_xyz)

        yy, xx, zz = sorted(s_xyz, key=lambda x: abs(x))

        aa = zz / 2.
        ar = (yy - xx) / 3.

        Saupe_components['Sxx'] = xx
        Saupe_components['Syy'] = yy
        Saupe_components['Szz'] = zz

        Saupe_components['Aa'].append(aa)
        Saupe_components['Ar'].append(ar)
        Saupe_components['Rh'].append(ar / abs(aa))

    # Calculate the summary statistics
    stats = calc_summary(magnetic_interactions, Saupe_components, data,
                         data_pred)

    return (data_pred, Saupe_components, stats)


# TODO: Implement DA running average to find dynamics (blocks of 10)