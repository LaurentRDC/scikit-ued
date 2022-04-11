# -*- coding: utf-8 -*-

from skued import autocenter, gaussian, kinematicsim
from crystals import Crystal
from scipy.ndimage import gaussian_filter
import numpy as np
import itertools as it
import pytest

np.random.seed(23)

DIFF_PATTERN_SIZE = 256


def diff_pattern_sc(center):
    """Simulate a single-crystal diffraction pattern"""
    r, c = center

    cryst = Crystal.from_database("BaTiO3_cubic")
    rr, cc = np.indices((DIFF_PATTERN_SIZE, DIFF_PATTERN_SIZE))
    rr -= r
    cc -= c

    kx = 10 * rr / (DIFF_PATTERN_SIZE // 2)
    ky = 10 * cc / (DIFF_PATTERN_SIZE // 2)

    kk = np.sqrt(kx**2 + ky**2)
    I = kinematicsim(crystal=cryst, kx=kx, ky=ky, energy=50)
    I[kk < 1] = 0
    gaussian_filter(I, sigma=2, output=I)
    I = np.minimum(I, 0.8 * I.max())
    I /= I.max()
    return I


def test_autocenter_trivial():
    """Test that autocenter() finds the center of perfect gaussian at (0,0)"""
    im = np.zeros(shape=(256, 256), dtype=float)
    center = np.asarray(im.shape) / 2
    rows, cols = np.indices(im.shape)
    im += gaussian([rows, cols], center=center, fwhm=center.max() / 2)

    assert np.allclose(autocenter(im), center, atol=1)


def test_autocenter_no_side_effects():
    """Test that autocenter() does not modify the inputs"""
    im = np.random.random(size=(256, 256))
    mask = np.ones_like(im, dtype=bool)

    # Modifying the arrays will result in "ValueError: output array is read-only"
    im.setflags(write=False)
    mask.setflags(write=False)

    autocenter(im, mask=mask)


@pytest.mark.parametrize("rc", range(-10, 10, 2))
@pytest.mark.parametrize("cc", range(-10, 10, 2))
def test_autocenter_gaussian_shifted(rc, cc):
    """Test that autocenter() finds the center of a shifted gaussian"""
    im = np.zeros(shape=(128, 128), dtype=float)
    rows, cols = np.indices(im.shape)
    center = np.array([64 + rc, 64 + cc])
    im += gaussian([rows, cols], center=center, fwhm=20)

    assert np.allclose(autocenter(im), center, atol=1)


@pytest.mark.parametrize("rc", range(-10, 10, 5))
@pytest.mark.parametrize("cc", range(-10, 10, 5))
def test_autocenter_shifted_with_mask(rc, cc):
    """Test that autocenter() finds the center of a shifted gaussian, where
    the center has been masked away."""
    im = np.zeros(shape=(256, 256), dtype=float)
    mask = np.ones_like(im, dtype=bool)
    mask[0:130, 118:138] = False

    rows, cols = np.indices(im.shape)
    center = np.array([128 + rc, 128 + cc])
    im += gaussian([rows, cols], center=center, fwhm=50)

    im[np.logical_not(mask)] *= 0.8

    assert np.allclose(autocenter(im, mask=mask), center, atol=1)


CENTERS = list(
    range(DIFF_PATTERN_SIZE // 3, 2 * DIFF_PATTERN_SIZE // 3, DIFF_PATTERN_SIZE // 10)
)


@pytest.mark.parametrize("rc", CENTERS)
@pytest.mark.parametrize("cc", CENTERS)
def test_autocenter_single_crystal(rc, cc):
    """Test that autocenter() finds the center of a simulated
    shifted single-crystal diffraction pattern."""
    I = 10 * diff_pattern_sc(center=(rc, cc))
    mask = np.ones_like(I, dtype=bool)
    mask[0 : rc + 10, cc - 10 : cc + 10] = False
    I += 0.01 * I.max() * np.random.random(size=I.shape)

    assert np.allclose(autocenter(I, mask=mask), (rc, cc), atol=1)


@pytest.mark.parametrize("rc", CENTERS)
@pytest.mark.parametrize("cc", CENTERS)
def test_autocenter_single_crystal_ewald_walkoff(rc, cc):
    """Test that autocenter() finds the center of a simulated
    shifted single-crystal diffraction pattern."""
    I = 10 * diff_pattern_sc(center=(rc, cc))

    # Walkoff is the effect when diffraction patterns slide reflections
    # at strange angles, such that certain bragg peaks that should have the same
    # intensity do not (e.g. brigher (n00) vs (-n00))
    # We simulate this with a linear intensity gradient
    rows, cols = np.indices(I.shape)
    walkoff = rows.astype(float) / rows.max()
    walkoff *= 0.2
    I += walkoff

    mask = np.ones_like(I, dtype=bool)
    mask[0 : rc + 10, cc - 10 : cc + 10] = False
    I += 0.01 * I.max() * np.random.random(size=I.shape)

    assert np.allclose(autocenter(I, mask=mask), (rc, cc), atol=1)
