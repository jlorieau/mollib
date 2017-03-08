from math import sin, cos, pi, sqrt

import numpy as np

def R_euler_zyz(alpha, beta, gamma):
    """Generate the Euler Z-Y-Z rotation matrix.
    """
    # Trig.
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

def R(axis, theta):
    """Generate the rotation matrix about the given axis.
    """
    axis = np.asarray(axis)
    axis = axis / sqrt(np.dot(axis, axis))
    a = cos(theta * pi / 360.0)
    b, c, d = -axis * sin(theta * pi / 360.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
