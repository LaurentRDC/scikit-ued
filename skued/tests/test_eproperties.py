# -*- coding: utf-8 -*-


import numpy as np

from skued import electron_wavelength, electron_velocity, interaction_parameter, lorentz
from scipy.constants import speed_of_light


def test_lorentz_trivial():
    """Test that the lorentz factor for 0 kV is unity"""
    assert lorentz(0) == 1


def test_lorentz_vectorized():
    """Test lorentz() on an array of energies"""
    kV = np.zeros((128,), dtype=float)
    assert np.allclose(lorentz(kV), np.ones_like(kV))


def test_lorentz_range():
    """Test that lorentz factor is always in the range [1, infty)"""
    kv = np.linspace(0, 1e6, num=256)
    factors = lorentz(kv)
    assert np.all(np.greater_equal(factors, 1))


def test_electron_wavelength_known():
    """Test that the electron wavelength at certain voltages is as expected."""
    assert round(abs(electron_wavelength(10) - 0.122), 3) == 0
    assert round(abs(electron_wavelength(200) - 0.0250), 3) == 0


def test_electron_velocity_trivial():
    """Test that the electron velocity at zero energy is zero"""
    assert electron_velocity(0) == 0


def test_electron_velocity_limits():
    """Test that the electron velocity never exceeds the speed of light."""
    c = speed_of_light * 1e10  # Speed of light in Ang/s
    assert electron_velocity(1e20) / c == 1


def test_interaction_parameter_100kV():
    """Test that the interaction_parameter(100) is what is expected
    from Kirkland 2010"""
    assert round(abs(interaction_parameter(100) - 0.924 * 1e-3), 6) == 0


def test_interaction_parameter_vectorized():
    """Test that interaction_parameters() is vectorized"""
    kV = np.full((64,), fill_value=100.0)
    i = np.full_like(kV, fill_value=0.924 * 1e-3)
    assert np.allclose(interaction_parameter(kV), i, atol=1e-6)
