# -*- coding: utf-8 -*-

import unittest
from .. import mad
import numpy as np

class TestMad(unittest.TestCase):

    def test_trivial(self):
        """ Test that the median absolute dev of an array of zeroes is zero """
        arr = np.zeros((16,))
        self.assertTrue(np.allclose(mad(arr), 0))
    
    def test_integers(self):
        """ Test that mad of an integer array is working as intended """
        arr_int = np.random.randint(-15, 15, size = (8,8))
        arr_flo = np.array(arr_int, copy = True, dtype = np.float)

        mad_int = mad(arr_int)
        mad_flo = mad(arr_flo)

        self.assertTrue(np.allclose(mad_int, mad_flo))
    
    def test_correctness(self):
        """ Test that the mad() function works as expected """
        arr = np.random.random(size = (32,))
        from_mad = mad(arr)
        explicit = np.abs(arr - np.median(arr))

        self.assertTrue(np.allclose(from_mad, explicit))
    
    def test_side_effects(self):
        """ Test that input array is not modified by mad() """
        arr = np.random.random(size = (32,))
        arr.setflags(write = False)
        mad(arr)

if __name__ == '__main__':
    unittest.main()