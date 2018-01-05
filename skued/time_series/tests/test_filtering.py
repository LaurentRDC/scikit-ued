# -*- coding: utf-8 -*-

import unittest
from .. import mad, outliers_mad
import numpy as np

class TestMad(unittest.TestCase):

    def test_trivial(self):
        """ Test that the median absolute dev of an array of zeroes is zero """
        arr = np.zeros((16,))
        self.assertEqual(mad(arr), 0)
    
    def test_correctness(self):
        """ Test that the mad() function works as expected """
        arr = np.random.random(size = (32,))
        from_mad = mad(arr)
        explicit = 1.4826 * np.median(np.abs(arr - np.median(arr)))

        self.assertTrue(np.allclose(from_mad, explicit))

class TestOutliersMad(unittest.TestCase):
    
    def test_trivial(self):
        """ Test that all elements of an array are outliers if the
        threshold is zero """
        arr = np.random.random(size = (64,))
        mask = outliers_mad(arr, thresh = 0)
        self.assertTrue(np.all(arr))


if __name__ == '__main__':
    unittest.main()