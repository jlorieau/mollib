# -*- coding: utf-8 -*-
"""
Functions to conduct the SVD on the RDCs and RACS.
"""
import numpy as np
from scipy import linalg

def Saupe_matrices(magnetic_interactions, data):
    """Calculate the best-fit Saupe matrices for the given interaction_arrays
    and data.

    Parameters
    ----------
    magnetic_interactions: list of dicts
        A list of dicts for each molecule to be fit.
        Each dict consists of (key) the interaction label (ex: '14N-H') and
        (value) a tuple of the scaling constant and the 5x1 array of the
        interaction for the dipolar or chemical shift/
        (Generated by the process_molecule Process objects)

    data: dicts
        A dict with the tensor key as keys, and a RDC
        or RACS datum as values.

    Returns
    -------
    (data_pred, S_xyz, Aa, Ar, Rh): tuple
        data_pred: dict
            A dict with the (key) interaction labels and (value) predicted RDCs
            and RACSs.
        S_xyz: The 3x1 Saupe matrix in x/y/z representation.
        Aa: The degree of alignment.
        Ar: The alignment rhombicity.
        Rh: The rhombicity

    """
    assert isinstance(magnetic_interactions, list)

    # Make an ordered list of the keys
    ordered_keys = sorted(data.keys())

    # Construct the A and D matrices
    A = []
    D = []
    for key in ordered_keys:
        for interaction_dict in magnetic_interactions:
            if key not in interaction_dict:
                msg = ("The orientation for data point {} has not been "
                       "calculated.")
                raise Exception(msg.format(key))
            scale, arr = interaction_dict[key]
            A.extend([arr,])

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
        data_pred[key] = D

    # Break up the S matrix into individual Saupe matrices, Das and Rh
    S_xyz = []
    Aa, Ar, Rh = [], [], []
    for x in range(0, len(S), 5):
        s = S[x:x + 5]
        s_xyz = np.array([[-0.5 * (s[0] - s[1]), s[2], s[3], ],
                         [s[2], -0.5 * (s[0] + s[1]), s[4]],
                         [s[3], s[4], s[0]]])
        s_xyz = linalg.eigvals(s_xyz).real
        S_xyz.append(s_xyz)

        xx, yy, zz = [i for i in sorted(abs(s_xyz))]
        da = max(s_xyz) / 2. if max(s_xyz) == zz else min(s_xyz) / 2.
        dr = (yy - xx) / 3.
        Aa.append(da)
        Ar.append(dr)
        Rh.append(dr / abs(da))

    return data_pred, S_xyz, Aa, Ar, Rh


