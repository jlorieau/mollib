# -*- coding: utf-8 -*-
"""
Functions to conduct the SVD on the RDCs and RACS.
"""
import numpy as np
from scipy import linalg

def Saupe_matrices(interaction_arrays, data):
    """Calculate the best-fit Saupe matrices for the given interaction_arrays
    and data.

    Parameters
    ----------
    interaction_arrays: list of dicts
        A list of dicts for each molecule to be fit.
        Each dict consists of the interaction label (ex: '14N-H') as keys and a
        5x1 array of the interaction values for the dipolar or chemical shift,
        generated by the process_molecule Process objects, as values.

    data: dicts
        A dict with the tensor key as keys, and a RDC
        or RACS datum as values.

    Returns
    -------
    Saupe_list: list
        A list of 5x1 Saupe matrices for each molecule/conformer.
    """
    assert isinstance(interaction_arrays, list)

    # Make an ordered list of the keys
    ordered_keys = sorted(data.keys())

    # Construct the A and D matrices
    A = []
    D = []
    for key in ordered_keys:
        for interaction_dict in interaction_arrays:
            if key not in interaction_dict:
                msg = ("The orientation for data point {} has not been "
                       "calculated.")
                raise Exception(msg.format(key))
            A.extend([interaction_dict[key]])

            expt_value = data[key].value
            expt_error = data[key].error
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

    # Calculate the Saupe matrix
    A_inv = np.dot(V.transpose(), np.dot(w_inv, U.transpose()))
    Saupe = np.dot(A_inv, D)

    # Calculate the A matrix for all the interactions in the interaction_dict.
    keys = {i for d in interaction_arrays for i in d.keys()}
    ordered_keys = sorted(keys)
    A = []
    for key in ordered_keys:
        for interaction_dict in interaction_arrays:
            if key not in interaction_dict:
                continue
            A.extend([interaction_dict[key]])

    D_pred = np.dot(A, Saupe)

    # Calculate the predicted data. It must be sorted the same way as the
    # interaction_dict
    data_pred = {}
    for key, D in zip(ordered_keys, D_pred):
        data_pred[key] = D

    return Saupe, data_pred

