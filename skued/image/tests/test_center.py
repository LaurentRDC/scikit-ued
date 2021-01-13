# -*- coding: utf-8 -*-

from skued import autocenter, gaussian
import numpy as np
import itertools as it
import pytest

np.random.seed(23)


def test_autocenter_trivial():
    """ Test that autocenter() finds the center of perfect gaussian at (0,0) """
    im = np.zeros(shape=(128, 128), dtype=np.float)
    center = np.asarray(im.shape) / 2
    rows, cols = np.indices(im.shape)
    im += gaussian([rows, cols], center=center, fwhm=center.max() / 2)

    assert np.allclose(autocenter(im), center)


@pytest.mark.parametrize("rc", range(-10, 10, 2))
@pytest.mark.parametrize("cc", range(-10, 10, 2))
def test_autocenter_shifted(rc, cc):
    """ Test that autocenter() finds the center of a shifted gaussian """
    im = np.zeros(shape=(128, 128), dtype=np.float)
    rows, cols = np.indices(im.shape)
    center = np.array([64 + rc, 64 + cc])
    im += gaussian([rows, cols], center=center, fwhm=20)

    assert np.allclose(autocenter(im), center, atol=1)


@pytest.mark.parametrize("rc", range(-10, 10, 5))
@pytest.mark.parametrize("cc", range(-10, 10, 5))
def test_autocenter_shifted_with_mask(rc, cc):
    """Test that autocenter() finds the center of a shifted gaussian, where
    the center has been masked away."""
    im = np.zeros(shape=(256, 256), dtype=np.float)
    mask = np.ones_like(im, dtype=np.bool)
    mask[0:130, 118:138] = False

    rows, cols = np.indices(im.shape)
    center = np.array([128 + rc, 128 + cc])
    im += gaussian([rows, cols], center=center, fwhm=50)

    im[np.logical_not(mask)] *= 0.8

    assert np.allclose(autocenter(im, mask=mask), center, atol=1)
