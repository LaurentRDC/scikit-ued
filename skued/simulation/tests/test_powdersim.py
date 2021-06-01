# -*- coding: utf-8 -*-

from copy import deepcopy

import numpy as np

from crystals import Crystal
from skued import powdersim


def test_powdersim_return_shape():
    """Test that the return shape of powdersim() is as expected"""
    q = np.linspace(2, 10, 200)
    pattern = powdersim(Crystal.from_database("C"), q)
    assert pattern.shape == q.shape


def test_powdersim_peak_alignment():
    """Test that the diffraction peaks align with what is expected."""
    crystal = Crystal.from_database("C")

    for reflection in [(0, 1, 1), (1, 2, 0), (-1, 2, 0)]:
        qknown = np.linalg.norm(crystal.scattering_vector((0, 1, 1)))

        # Range of scattering vectors is tightly centered around a particular reflection
        # So that the maximum of the diffraction pattern MUST be at reflection (010)
        q = np.linspace(qknown - 0.1, qknown + 0.1, 256)
        pattern = powdersim(Crystal.from_database("C"), q)
        assert abs(q[np.argmax(pattern)] - qknown) < q[1] - q[0]
