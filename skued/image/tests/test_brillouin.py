# -*- coding: utf-8 -*-

from skued import kinematicsim, brillouin_zones, autocenter, bragg_peaks_persistence
from crystals import Crystal
from scipy.ndimage import gaussian_filter
import numpy as np
import itertools as it
import pytest

DIFF_PATTERN_SIZE = 256
MAX_INV_ANG = 5.5


def diff_pattern_sc():
    """Simulate a single-crystal diffraction pattern"""
    np.random.seed(23)
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


def test_brillouin_consistency():
    kx, ky, I, cryst = diff_pattern_sc()
    peaks, _, _, _ = bragg_peaks_persistence(I, prominence=0.04)
    center = autocenter(I)
    BZ = brillouin_zones(
        I,
        mask=np.ones_like(I, dtype=bool),
        peaks=peaks.astype(int),
        center=center.astype(int),
    )
    BZ.getVisibleBZs()
    assert 1 - BZ.determineConsistency() < 0.03
