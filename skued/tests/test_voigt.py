# -*- coding: utf-8 -*-
import numpy as np
import random
import pytest

from skued import gaussian, lorentzian, pseudo_voigt
from skued.voigt import _pseudo_voigt_mixing_factor


def integrate_1d(x, f):
    """Numerically integrate a function f(x)."""
    return np.trapz(f, x)


def integrate_2d(x, y, f):
    """Numerically-integrate a function f(x, y)."""
    return np.trapz(np.trapz(f, y[None, :], axis=1), x, axis=0)


def integrate_3d(x, y, z, f):
    """Numerically-integrate a function f(x, y)."""
    return np.trapz(
        np.trapz(np.trapz(f, z[None, None, :], axis=2), y[None, :], axis=1), x, axis=0
    )


multifunc = pytest.mark.parametrize(
    "func", (gaussian, lorentzian, lambda x, c, w: pseudo_voigt(x, c, w, w))
)


@multifunc
def test_1d_unit_integral(func):
    dx = 0.1
    x = np.arange(-10, 10, dx)
    center = 0

    integral = integrate_1d(x=x, f=func(x, center, 0.2))
    assert round(abs(integral - 1), 2) == 0


@multifunc
def test_2d_unit_integral(func):
    dx = dy = 0.1
    extent = np.arange(-10, 10, dx)
    xx, yy = np.meshgrid(extent, extent)
    center = [0, 0]

    integral = integrate_2d(extent, extent, func([xx, yy], center, 0.2))
    assert round(abs(integral - 1), 2) == 0


@multifunc
def test_3d_unit_integral(func):
    dx = dy = dz = 0.01
    extent = np.arange(-1.3, 1.3, dx)
    xx, yy, zz = np.meshgrid(extent, extent, extent)
    center = [0, 0, 0]

    integral = integrate_3d(extent, extent, extent, func([xx, yy, zz], center, 0.1))
    assert round(abs(integral - 1), 2) == 0


def test_bounds_pseudo_mixing_factor():
    """Test that the mixing factor is always between 0 and 1"""
    random.seed(0)
    for _ in range(100):
        fl = random.uniform(1e-8, 100)
        fg = random.uniform(1e-8, 100)
        eta = _pseudo_voigt_mixing_factor(fl, fg)
        assert eta >= 0
        assert eta <= 1
