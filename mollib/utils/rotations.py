"""
Utility function for generating rotation matrices, getting rotation angles and
dealing with rotations.
"""

from math import sin, cos, atan2, pi, sqrt

import numpy as np


def R_euler_zyz(alpha, beta, gamma):
    """Generate the Euler Z-Y-Z rotation matrix.

    .. note::

        This is a right-handed rotation of the frame::

              z |
                |
                |
                |
                -------- y
               /
              /
           x /


    Parameters
    ----------
    alpha: float
        The first rotation about the z-axis (in degrees)
    beta: float
        The second rotation about the new y-axis (in degrees)
    gamma: float
        The third rotation about the new z-axis (in degrees)

    Returns
    -------
    rotation_matrix: array
        The 3x3 rotation matrix. (Row vectors)
    """
    sin_a = sin(alpha * pi / 180.)
    sin_b = sin(beta * pi / 180.)
    sin_g = sin(gamma * pi / 180.)

    cos_a = cos(alpha * pi / 180.)
    cos_b = cos(beta * pi / 180.)
    cos_g = cos(gamma * pi / 180.)

    return [[-sin_a * sin_g + cos_a * cos_b * cos_g,
             cos_a * sin_g + sin_a * cos_b * cos_g,
             -sin_b * cos_g],
            [-sin_a * cos_g - cos_a * cos_b * sin_g,
             cos_a * cos_g - sin_a * cos_b * sin_g,
             sin_b * sin_g],
            [cos_a * sin_b,
             sin_a * sin_b,
             cos_b]]


def euler_zyz(R):
    """Generate the Z-Y-Z Euler angles (in degrees) for the given rotation
    matrix.

    Parameters
    ----------
    R: array
        The 3x3 rotation matrix.

    Returns
    -------
    angles: tuple
        The alpha, beta and gamma angles of the rotation.

        - 'alpha': The first rotation angle, about the z-axis (in degrees).
          [0, 360]
        - 'beta': The second rotation angle, about the new y-axis (in degrees).
          [0, 180]
        - 'gamma: The third rotation angle, about the new z-axis (in degrees).
          [0, 360]

        .. note:: The Euler angles are subject to a gimbal lock in which the
                  alpha and gamma angles become interdependent when beta is
                  near 0 degrees . In this case, the azimuthal angle will be
                  reported to alpha (and gamma=0 degrees).

    Examples
    --------
    >>> import numpy as np
    >>> m = R_euler_zyz(0., 90, 0)
    >>> angles = np.round(euler_zyz(m), 1)
    >>> print(tuple(angles))
    (0.0, 90.0, 0.0)
    >>> m = R_euler_zyz(0., 10, 90.)
    >>> angles = np.round(euler_zyz(m), 1)
    >>> print(tuple(angles))
    (0.0, 10.0, 90.0)
    >>> m = R_euler_zyz(90., 10, 0.)
    >>> angles = np.round(euler_zyz(m), 1)
    >>> print(tuple(angles))
    (90.0, 10.0, 0.0)
    """
    sb = sqrt(R[2][0] ** 2 + R[2][1] ** 2)

    beta = atan2(sb, R[2][2]) * 180. / pi
    if sb < 0.001:
        # This is the gimbol lock situation. When beta is very small, then
        # alpha and gamma rotations are effectively coupled, and only one can
        # be defined.
        alpha = atan2(R[0][1], R[0][0]) * 180. / pi
        gamma = 0.0
    else:
        alpha = atan2(R[2][1], R[2][0]) * 180. / pi
        gamma = atan2(R[1][2], -R[0][2]) * 180. / pi

    # Return angles between [0, 360] or [0, 180]
    alpha = (alpha + 360) % 360
    beta = (beta + 180) % 180
    gamma = (gamma + 360) % 360

    return alpha, beta, gamma


def R(axis, theta):
    """Generate the rotation matrix about the given axis.

    .. note::

        This is a right-handed rotation of the frame::

              z |
                |
                |
                |
                -------- y
               /
              /
           x /


    Parameters
    ----------
    axis: array
        The 3x1 vector to conduct the rotation about. (This vector will be
        normalized.)
    theta: float
        The rotation angle (in degrees).

    Returns
    -------
    rotation_matrix: array
        The 3x3 rotation matrix.

    Examples
    --------
    >>> import numpy as np
    >>> Rz = lambda theta: R([0., 0., 1], theta)
    >>> list(map(lambda x: round(x, 1), np.dot(Rz(90.), [1., 0., 0.])))
    [0.0, 1.0, 0.0]
    >>> list(map(lambda x: round(x, 1), np.dot(Rz(90.), [0., 1., 0.])))
    [-1.0, 0.0, 0.0]
    >>> list(map(lambda x: round(x, 1), np.dot(Rz(90.), [0., 0., 1.])))
    [0.0, 0.0, 1.0]
    """
    axis = np.asarray(axis)
    axis /= sqrt(np.dot(axis, axis))
    a = cos(theta * pi / 360.0)
    b, c, d = -axis * sin(theta * pi / 360.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

