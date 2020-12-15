# -*- coding: utf-8 -*-
import unittest

import numpy as np
import random

from skued import gaussian, lorentzian, pseudo_voigt
from skued.voigt import _pseudo_voigt_mixing_factor


def integrate_1d(x, f):
    """ Numerically integrate a function f(x)."""
    return np.trapz(f, x)


def integrate_2d(x, y, f):
    """ Numerically-integrate a function f(x, y). """
    return np.trapz(np.trapz(f, y[None, :], axis=1), x, axis=0)


def integrate_3d(x, y, z, f):
    """ Numerically-integrate a function f(x, y). """
    return np.trapz(
        np.trapz(np.trapz(f, z[None, None, :], axis=2), y[None, :], axis=1), x, axis=0
    )


class UnitIntegral(object):
    # Same tests for lorentzian and gaussian
    # Unit integral testing is only to one decimal place
    # to save some time.

    def test_1d_unit_integral(self):
        dx = 0.1
        x = np.arange(-10, 10, dx)
        center = 0

        integral = integrate_1d(x=x, f=self.function(x, center, 0.2))
        self.assertAlmostEqual(integral, 1, places=2)

    def test_2d_unit_integral(self):
        dx = dy = 0.1
        extent = np.arange(-10, 10, dx)
        xx, yy = np.meshgrid(extent, extent)
        center = [0, 0]

        integral = integrate_2d(extent, extent, self.function([xx, yy], center, 0.2))
        self.assertAlmostEqual(integral, 1, places=2)

    def test_3d_unit_integral(self):
        dx = dy = dz = 0.01
        extent = np.arange(-1.3, 1.3, dx)
        xx, yy, zz = np.meshgrid(extent, extent, extent)
        center = [0, 0, 0]

        integral = integrate_3d(
            extent, extent, extent, self.function([xx, yy, zz], center, 0.1)
        )
        self.assertAlmostEqual(integral, 1, places=2)


class TestGaussianUnitIntegral(UnitIntegral, unittest.TestCase):
    def setUp(self):
        self.function = gaussian


class TestLorentzianUnitIntegral(UnitIntegral, unittest.TestCase):
    def setUp(self):
        self.function = lorentzian


class TestVoigtFunctionUnitIntegral(UnitIntegral, unittest.TestCase):
    # voigt and pseudo_voigt have two width parameters
    # Unit integral testing is only to one decimal place
    # to save some time.

    def setUp(self):
        self.function = lambda x, c, w: pseudo_voigt(x, c, w, w)


class TestPseudoMixingFactor(unittest.TestCase):
    """ Test the pseudo-mixing factor "eta" that determines the fraction of Gaussian vs Lorentzian. """

    def test_bounds(self):
        """ Test that the mixing factor is always between 0 and 1 """
        random.seed(0)
        for _ in range(100):
            fl = random.uniform(1e-8, 100)
            fg = random.uniform(1e-8, 100)
            eta = _pseudo_voigt_mixing_factor(fl, fg)
            self.assertGreaterEqual(eta, 0)
            self.assertLessEqual(eta, 1)


if __name__ == "__main__":
    unittest.main()
