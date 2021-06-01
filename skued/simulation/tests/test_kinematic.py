# -*- coding: utf-8 -*-


import numpy as np

from skued.simulation.kinematic import fft2freq, limit_bandwidth
import pytest


def test_shape():
    """Test that the output from fft2freq has the same shape as the input."""
    extent_x = np.linspace(0, 1)
    extent_y = np.linspace(3, 4, num=47)

    for indexing in {"ij", "xy"}:
        x, y = np.meshgrid(extent_x, extent_y, indexing=indexing)
        kx, ky = fft2freq(x, y, indexing=indexing)
        assert x.shape == kx.shape
        assert y.shape == ky.shape


def test_indexing_error():
    """Test that fft2freq correctly raises ValueError for invalid indexing."""
    extent_x = np.linspace(0, 1)
    extent_y = np.linspace(3, 4, num=47)

    with pytest.raises(ValueError):
        fft2freq([1, 2], [1, 2], indexing="ab")


def test_vs_fftfreq():
    """Test that the results make sense with respect to 1D case"""
    extent_x = np.arange(0, 10, step=0.1)
    extent_y = np.arange(3, 4, step=0.1)

    for indexing in {"ij", "xy"}:
        x, y = np.meshgrid(extent_x, extent_y, indexing=indexing)
        kx, ky = fft2freq(x, y, indexing=indexing)

        # same array creates from the 1d case
        kx_1d = np.fft.fftfreq(len(extent_x), d=0.1)
        ky_1d = np.fft.fftfreq(len(extent_y), d=0.1)
        kx2, ky2 = np.meshgrid(kx_1d, ky_1d, indexing=indexing)

        assert np.allclose(kx, kx2)
        assert np.allclose(ky, ky2)


def test_idempotence():
    """Test that applying limit_bandwidth more than once has no effect"""
    im = np.ones((32, 32), dtype=float)
    kx, ky = np.meshgrid(np.arange(-16, 16), np.arange(-16, 16))
    k = np.hypot(kx, ky)

    limited1 = limit_bandwidth(im, k, k.max() / 2)
    limited2 = limit_bandwidth(limited1, k, k.max() / 2)

    assert np.allclose(limited1, limited2)
