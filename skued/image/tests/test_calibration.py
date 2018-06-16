# -*- coding: utf-8 -*-

import numpy as np
from .. import powder_calq
from ... import Crystal, powdersim
import unittest

from skimage.filters import gaussian

np.random.seed(23)

class TestPowderCalQ(unittest.TestCase):

    def test_simulation(self):
        """ 
        Test calibration from simulation, down to 1% error . Peaks (200) and (220) from monoclinic VO2 are used,
        """

        s = np.linspace(0.11, 0.8, 1024)
        q = 4*np.pi*s
        c = Crystal.from_database('vo2-m1')
        I = powdersim(c, s)

        peak1 = (2, 0, 0)
        Gx1, Gy1, Gz1 = c.scattering_vector(*peak1)
        q1 = np.sqrt(Gx1**2 + Gy1**2 + Gz1**2)
        arr_index1 = np.argmin(np.abs(q - q1))

        peak2 = (2, 2, 0)
        Gx2, Gy2, Gz2 = c.scattering_vector(*peak2)
        q2 = np.sqrt(Gx2**2 + Gy2**2 + Gz2**2)
        arr_index2 = np.argmin(np.abs(q - q2))

        calibrated = powder_calq(I, c, peak_indices = (arr_index1, arr_index2), miller_indices = (peak1, peak2))

        self.assertTupleEqual(I.shape, calibrated.shape)
        self.assertTrue(np.allclose(q, calibrated, rtol = 0.01))

    def test_simulation_3_peaks(self):
        """ 
        Test calibration from simulation, down to 1% error . Peaks (200) and (220) from monoclinic VO2 are used,
        """
        s = np.linspace(0.11, 0.8, 1024)
        q = 4*np.pi*s
        c = Crystal.from_database('vo2-m1')
        I = powdersim(c, s)

        peak1 = (2, 0, 0)
        Gx1, Gy1, Gz1 = c.scattering_vector(*peak1)
        q1 = np.sqrt(Gx1**2 + Gy1**2 + Gz1**2)
        arr_index1 = np.argmin(np.abs(q - q1))

        peak2 = (2, 2, 0)
        Gx2, Gy2, Gz2 = c.scattering_vector(*peak2)
        q2 = np.sqrt(Gx2**2 + Gy2**2 + Gz2**2)
        arr_index2 = np.argmin(np.abs(q - q2))

        peak3 = (3, 0, -2)
        Gx2, Gy2, Gz2 = c.scattering_vector(*peak3)
        q3 = np.sqrt(Gx2**2 + Gy2**2 + Gz2**2)
        arr_index3 = np.argmin(np.abs(q - q3))

        calibrated = powder_calq(I, c, peak_indices = (arr_index1, arr_index2, arr_index3), miller_indices = (peak1, peak2, peak3))

        self.assertTupleEqual(I.shape, calibrated.shape)
        self.assertTrue(np.allclose(q, calibrated, rtol = 0.01))

if __name__ == '__main__':
    unittest.main()