# -*- coding: utf-8 -*-
import unittest
from copy import deepcopy

import numpy as np

from crystals import Crystal
from skued import powdersim


class TestPowdersim(unittest.TestCase):
    def test_return_shape(self):
        """ Test that the return shape of powdersim() is as expected """
        q = np.linspace(2, 10, 200)
        pattern = powdersim(Crystal.from_database("C"), q)
        self.assertSequenceEqual(pattern.shape, q.shape)

    def test_peak_alignment(self):
        """ Test that the diffraction peaks align with what is expected. """
        crystal = Crystal.from_database("C")

        for reflection in [(0, 1, 1), (1, 2, 0), (-1, 2, 0)]:
            with self.subTest(f"Reflection {reflection}"):
                qknown = np.linalg.norm(crystal.scattering_vector((0, 1, 1)))

                # Range of scattering vectors is tightly centered around a particular reflection
                # So that the maximum of the diffraction pattern MUST be at reflection (010)
                q = np.linspace(qknown - 0.1, qknown + 0.1, 256)
                pattern = powdersim(Crystal.from_database("C"), q)
                self.assertAlmostEqual(q[np.argmax(pattern)], qknown, delta=q[1] - q[0])


if __name__ == "__main__":
    unittest.main()
