"""
Unit tests for the rotation utility functions.
"""
import unittest

import numpy as np

from mollib.utils.rotations import R_euler_zyz, euler_zyz

class TestRotations(unittest.TestCase):

    def test_R_euler_zyz(self):
        """Test the R_euler_zyz function."""

        # Start with no rotation
        m1 = np.round(R_euler_zyz(0., 0., 0.), 3)
        m2 = np.array([[1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0],
                       [0.0, 0.0, 1.0]])
        self.assertTrue(np.array_equal(m1, m2))

        # Then conduct an individual rotation about every axis.
        m1 = np.round(R_euler_zyz(90., 0., 0.), 3)
        m2 = np.array([[0.0,  1.0, 0.0],   # x -> y
                       [-1.0, 0.0, 0.0],   # y -> -x
                       [0.0,  0.0, 1.0]])  # z -> z
        self.assertTrue(np.array_equal(m1, m2))

        m1 = np.round(R_euler_zyz(0., 90., 0.), 3)
        m2 = np.array([[0.0, 0.0, -1.0],   # x -> -z
                       [0.0, 1.0,  0.0],   # y -> y
                       [1.0, 0.0,  0.0]])  # z -> x
        self.assertTrue(np.array_equal(m1, m2))

        m1 = np.round(R_euler_zyz(0., 0., 90.), 3)
        m2 = np.array([[0.0,  1.0, 0.0],   # x -> y
                       [-1.0, 0.0, 0.0],   # y -> -x
                       [0.0,  0.0, 1.0]])  # z -> z
        self.assertTrue(np.array_equal(m1, m2))

        # Test two concurrent rotations
        m1 = np.round(R_euler_zyz(90., 0., -90.), 3)
        m2 = np.array([[1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0],
                       [0.0, 0.0, 1.0]])
        self.assertTrue(np.array_equal(m1, m2))

        m1 = np.round(R_euler_zyz(90., 90., 0.), 3)
        m2 = np.array([[0.0,  0.0, -1.0],   # x -> -z
                       [-1.0, 0.0,  0.0],   # y -> -x
                       [0.0,  1.0,  0.0]])  # z -> y
        self.assertTrue(np.array_equal(m1, m2))

        m1 = np.round(R_euler_zyz(0., 90., 90.), 3)
        m2 = np.array([[0.0, 1.0, 0.0],   # x -> y
                       [0.0, 0.0, 1.0],   # y -> z
                       [1.0, 0.0, 0.0]])  # z -> x
        self.assertTrue(np.array_equal(m1, m2))

    def test_euler_zyz(self):
        """Test the euler_zyz function."""

        def angles(alpha, beta, gamma):
            """Return the angles from a rotation matrix constructed from
            the given angles."""
            R = R_euler_zyz(alpha, beta, gamma)
            return euler_zyz(R)

        # Try for beta of 0.1-179.9 degrees
        beta_range = [0.1,] + [float(b) for b in range(30, 170, 30)] + [179.9,]
        for beta in beta_range:
            for alpha in range(0, 720, 30):
                for gamma in range(0, 720, 30):
                    alpha = float(alpha)
                    beta = float(beta)
                    gamma = float(gamma)
                    a, b, g = angles(alpha, beta, gamma)

                    self.assertAlmostEqual(a, (alpha + 360) % 360)
                    self.assertAlmostEqual(b, beta)
                    self.assertAlmostEqual(g, (gamma + 360) % 360)

        # Due to the gimbal lock, the alpha/gamma angles become degenerate
        # when beta ~= 0 or 180 degrees
        for alpha, beta, gamma in ((45., 0.0, 45.),
                                   (0., 0.01, 130.),
                                   ):
            a, b, g = angles(alpha, beta, gamma)

            self.assertAlmostEqual(a, alpha + gamma, 1)
            self.assertAlmostEqual(b, beta, 1)


