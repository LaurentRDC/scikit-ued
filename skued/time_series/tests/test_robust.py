# -*- coding: utf-8 -*-

import unittest
from .. import mad
import numpy as np

class TestMad(unittest.TestCase):

    def test_trivial(self):
        """ Test that the median absolute dev of an array of zeroes is zero """
        arr = np.zeros((16,))
        self.assertTrue(np.allclose(mad(arr), 0))
    
    def test_correctness(self):
        """ Test that the mad() function works as expected """
        arr = np.random.random(size = (32,))
        from_mad = mad(arr)
        explicit = np.abs(arr - np.median(arr))

        self.assertTrue(np.allclose(from_mad, explicit))

if __name__ == '__main__':
    unittest.main()