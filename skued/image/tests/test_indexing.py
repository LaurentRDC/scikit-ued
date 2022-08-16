# -*- coding: utf-8 -*-

from skued import kinematicsim, bragg_peaks
from crystals import Crystal
from scipy.ndimage import gaussian_filter
import numpy as np
import itertools as it
import pytest

DIFF_PATTERN_SIZE = 256
MAX_INV_ANG = 5.5


def diff_pattern_sc():
    """Simulate a single-crystal diffraction pattern"""
    cryst = Crystal.from_database("C")
    rr, cc = np.indices((DIFF_PATTERN_SIZE, DIFF_PATTERN_SIZE))
    rr -= rr.shape[1] // 2
    cc -= cc.shape[0] // 2

    kx = 5.5 * rr / (DIFF_PATTERN_SIZE // 2)
    ky = 5.5 * cc / (DIFF_PATTERN_SIZE // 2)

    kk = np.sqrt(kx**2 + ky**2)
    I = kinematicsim(crystal=cryst, kx=kx, ky=ky, energy=50)
    I[kk < 1] = 0
    gaussian_filter(I, sigma=1, output=I)
    I = np.minimum(I, 0.8 * I.max())
    I /= I.max()
    I += 0.05 * np.random.random(size=I.shape)
    return kx, ky, I, cryst


def test_bragg_peaks():
    """Test that the `bragg_peaks` function finds all Bragg peaks."""
    kx, ky, I, cryst = diff_pattern_sc()
    kk = np.sqrt(kx**2 + ky**2)

    # only in-plane refls (hk0) and not (000)
    # Also, some reflections will appear at the edge of the frame
    in_plane_refls = [
        refl
        for refl in cryst.bounded_reflections(kk.max())
        if (refl[2] == 0 and np.abs(cryst.scattering_vector(refl)[0]) < kx.max())
    ]

    peaks = bragg_peaks(I, mask=np.ones_like(I, dtype=bool))

    assert len(peaks) == len(in_plane_refls)


def test_bragg_peaks_no_mask():
    """Test that the `bragg_peaks` function finds all Bragg peaks, without supplying a mask file."""
    kx, ky, I, cryst = diff_pattern_sc()
    kk = np.sqrt(kx**2 + ky**2)

    # only in-plane refls (hk0) and not (000)
    # Also, some reflections will appear at the edge of the frame
    in_plane_refls = [
        refl
        for refl in cryst.bounded_reflections(kk.max())
        if (refl[2] == 0 and np.abs(cryst.scattering_vector(refl)[0]) < kx.max())
    ]

    peaks = bragg_peaks(I)

    assert len(peaks) == len(in_plane_refls)
