# -*- coding: utf-8 -*-
from ..wavelets import (dualtree, idualtree, dt_max_level, dualtree_wavelet, 
						dualtree_first_stage, 
						ALL_QSHIFT, ALL_FIRST_STAGE, ALL_COMPLEX_WAV)

import numpy as np
import pywt
import unittest

np.random.seed(23)


##############################################################################
###             COMPLEX WAVELET 
##############################################################################

class TestComplexWavelets(unittest.TestCase):

    def setUp(self):
        self.array = np.sin(np.arange(0, 10, step = 0.01))
    
    def test_first_stage(self):
        """ Test of perfect reconstruction of first stage wavelets. """
        for wavelet in ALL_FIRST_STAGE:
            for wav in dualtree_first_stage(wavelet):
                # Using waverec and wavedec instead of dwt and idwt because parameters
                # don't need as much parsing.
                self.assertTrue(np.allclose( self.array, pywt.waverec(pywt.wavedec(self.array, wav), wav) ))

##############################################################################
###           DUAL-TREE COMPLEX WAVELET TRANSFORM
##############################################################################

class TestDualTree(object):
    """ Skeleton for 1D and 2D testing. Tests are run from subclasses. """
    
    def test_perfect_reconstruction_level_0(self):
        coeffs = dualtree(data = self.array, level = 0)
        reconstructed = idualtree(coeffs = coeffs)
        self.assertTrue(np.allclose(self.array, reconstructed))
    
    def test_perfect_reconstruction_level_1(self):
        for first_stage in ALL_FIRST_STAGE:
            coeffs = dualtree(data = self.array, level = 1, first_stage = first_stage)
            reconstructed = idualtree(coeffs = coeffs, first_stage = first_stage)
            self.assertTrue(np.allclose(self.array, reconstructed))
    
    def test_perfect_reconstruction_multilevel(self):
        for first_stage in ALL_FIRST_STAGE:
            for wavelet in ALL_COMPLEX_WAV:
                for level in range(1, dt_max_level(data = self.array, first_stage = first_stage, wavelet = wavelet)):
                    coeffs = dualtree(data = self.array, level = level, first_stage = first_stage, wavelet = wavelet)
                    reconstructed = idualtree(coeffs = coeffs, first_stage = first_stage, wavelet = wavelet)
                    self.assertTrue(np.allclose(self.array, reconstructed))
    
    def test_axis(self):
        for axis in range(0, self.array.ndim):
            coeffs = dualtree(data = self.array, level = 2, axis = axis)
            reconstructed = idualtree(coeffs = coeffs, axis = axis)
            self.assertTrue(np.allclose(self.array, reconstructed))
    
    def test_axis_limits(self):
        with self.assertRaises(ValueError):
            coeffs = dualtree(data = self.array, level = 1, axis = self.array.ndim)

# Actual tests

class Test1D(TestDualTree, unittest.TestCase):
    def setUp(self):
        self.array = np.random.random(size = (100,))

class Test2D(TestDualTree, unittest.TestCase):
    def setUp(self):
        self.array = np.random.random(size = (50,50))

class Test3D(TestDualTree, unittest.TestCase):
    def setUp(self):
        self.array = np.random.random(size = (10,10,10))
    
if __name__ == '__main__':
    unittest.main()