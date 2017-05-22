# -*- coding: utf-8 -*-
"""
Functions to conduct the SVD on the RDCs and RACS.
"""
import logging

import numpy as np
from scipy import linalg

from mollib.utils.rotations import euler_zyz
from mollib.utils.interactions import sort_func
from . import logs
from .analysis import calc_summary
from .utils import get_data_type
from . import settings


def get_error(label, data):
    """Return the error for the given interaction.

    Parameters
    ----------
    label: str
        interaction label
    data: dict
        The experimental/observed RDC and RACS data.
        
        - **key**: interaction labels, str
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

    # Otherwise fetch a default value
    interaction_type = sort_func(label)[0]
    if interaction_type in settings.default_error:
        return settings.default_error[interaction_type]

    interaction_type_rev = '-'.join(interaction_type.split('-')[::-1])
    if interaction_type_rev in settings.default_error:
        return settings.default_error[interaction_type_rev]

    # Finally, see if there's a default value for the bond type.
    # The bond_type converts 'CA-CB' into 'C-C'
    bond_type = '-'.join([i[0] for i in interaction_type.split('-')])
    if bond_type in settings.default_error:
        return settings.default_error[bond_type]

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
        
        - **key**: interaction labels, str
        - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data values.

    Returns
    -------
    (data_pred, Saupe_components, stats): tuple
    
        - data_pred: dict
    
          - The SVD predicted RDC and RACS data.
          - **key**: interaction labels (str)
          - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data
            values.
            
        - Saupe_components: dict
        
          - 'S_xyz': The 3x1 Saupe matrix in x/y/z repr, list of arrays
          - 'Aa': The degree of alignment, list
          - 'Ar': The alignment rhombicity, list
          - 'Rh': The rhombicity, list
          - 'alpha_z': The ZYZ alpha angle (deg) of the Saupe matrix, float
          - 'beta_y': The ZYZ beta angle (deg) of the Saupe matrix, float
          - 'gamma_z': The ZYZ gamma angle (degrees) of the Saupe matrix, float
             
        - stats: dict
        
          - See :func:`calc_statistics`
    """
    assert(isinstance(magnetic_interactions, list))

    # Make an ordered list of the keys
    ordered_keys = sorted(data.keys())

    # Construct the A and D matrices
    A = []
    D = []

    for key in ordered_keys:
        # Get the experimental value and error
        expt_value = data[key].value
        expt_error = get_error(key, data)

        # Construct the A-matrix
        A_line = []
        for interaction_dict in magnetic_interactions:
            # Check to see if the interaction has been processed.
            if key not in interaction_dict.keys():
                if key not in logs.errors:
                        msg = ("Processing of data point '{}' is not "
                               "implemented.")
                        logging.error(msg.format(key))
                        logs.errors.add(key)
                continue

            scale, arr = interaction_dict[key]
            A_line.extend(arr * scale / expt_error)
        if A_line:
            A.append(A_line)
            D.append(expt_value / expt_error)

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

    # Calculate the predicted RDCs and RACSs. This is done by calculating
    # a new A-matrix, including interactions not in the data.
    keys_pred = {i for d in magnetic_interactions for i in d.keys()}
    data_pred = {}
    A_pred = []

    for key in keys_pred:
        expt_error = get_error(key, data)
        A_line = []
        for interaction_dict in magnetic_interactions:
            # Check to see if the interaction has been processed.
            if key not in interaction_dict:
                if key not in logs.errors:
                    msg = ("Processing of data point '{}' is not "
                           "implemented.")
                    msg = msg.format(key)
                    logging.error(msg)
                    logs.errors.add(key)
                continue

            scale, arr = interaction_dict[key]
            A_line.extend(arr * scale / expt_error)

        A_pred.append(A_line)

    D_pred = np.dot(A_pred, S)

    # Copy the predicted RDC and RACS values to the data_pred dict.
    for key, D in zip(keys_pred, D_pred):
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

    # Conduct the S-matrix calculations for each molecule:
    Saupe_components = {}
    for i in molecule_number:
        # Contruct the molecular identifier. This is empty ('') if there's one
        # molecule, or it's (1), (2), ... if there are multiple molecules
        id_ = ' ({})'.format(i) if i is not None else ''

        if i is None:
            s = S
        else:
            s = S[(i-1) * 5: (i-1) * 5 + 5]

        # The order of the Saupe matrix components are:
        # Cyy, Czz, Cxy, Cxz, Cyz
        s_xyz = np.array([[-s[0] - s[1], s[2], s[3]],
                          [s[2],         s[0], s[4]],
                          [s[3],         s[4], s[1]]])

        # Conduct the eigen value decomposition
        eig_values, eig_vectors = linalg.eig(s_xyz)

        # Get the Saupe matrix angles
        alpha, beta, gamma = euler_zyz(eig_vectors)
        Saupe_components['alpha_z' + id_] = alpha
        Saupe_components['beta_y' + id_] = beta
        Saupe_components['gamma_z' + id_] = gamma

        # Get the Saupe matrix values
        s_xyz = eig_values.real
        Saupe_components['S_xyz' + id_] = s_xyz

        yy, xx, zz = sorted(s_xyz, key=lambda x: abs(x))

        aa = zz / 2.
        ar = (yy - xx) / 3.

        Saupe_components['Sxx' + id_] = xx
        Saupe_components['Syy' + id_] = yy
        Saupe_components['Szz' + id_] = zz

        Saupe_components['Aa' + id_] = aa
        Saupe_components['Ar' + id_] = ar
        Saupe_components['Rh' + id_] = ar / aa

    # Calculate the summary statistics
    stats = calc_summary(magnetic_interactions, Saupe_components, data,
                         data_pred)

    return (data_pred, Saupe_components, stats)


# TODO: Implement DA running average to find dynamics (blocks of 10)