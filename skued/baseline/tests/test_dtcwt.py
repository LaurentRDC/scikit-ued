# -*- coding: utf-8 -*-
from ..dtcwt import (dtcwt, idtcwt, dt_max_level, dualtree_wavelet, 
						dt_first_stage, available_first_stage_filters, available_dt_filters)

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
        for wavelet in available_first_stage_filters():
            for wav in dt_first_stage(wavelet):
                # Using waverec and wavedec instead of dwt and idwt because parameters
                # don't need as much parsing.
                self.assertTrue(np.allclose( self.array, pywt.waverec(pywt.wavedec(self.array, wav), wav) ))

##############################################################################
###           DUAL-TREE COMPLEX WAVELET TRANSFORM
##############################################################################

class TestDtcwt(object):
    """ Skeleton for 1D and 2D testing. Tests are run from subclasses. """
    
    def test_perfect_reconstruction_level_0(self):
        """ Test perfect reconstruction for a 0-level decomposition """
        coeffs = dtcwt(data = self.array, level = 0, first_stage = 'sym6', wavelet = 'qshift1')
        reconstructed = idtcwt(coeffs = coeffs, first_stage = 'sym6', wavelet = 'qshift1')
        self.assertTrue(np.allclose(self.array, reconstructed))
    
    def test_perfect_reconstruction_level_1(self):
        """ Test perfect reconstruction for a single decomposition level """
        for first_stage in available_first_stage_filters():
            coeffs = dtcwt(data = self.array, level = 1, first_stage = first_stage, wavelet = 'qshift1')
            reconstructed = idtcwt(coeffs = coeffs, first_stage = first_stage, wavelet = 'qshift1')
            self.assertTrue(np.allclose(self.array, reconstructed))
    
    def test_perfect_reconstruction_multilevel(self):
        """ Test perfect reconstruction for all levels, for all first_stage wavelets, for all DT wavelets """
        for first_stage in available_first_stage_filters():
            for wavelet in available_dt_filters():
                for level in range(1, dt_max_level(data = self.array, first_stage = first_stage, wavelet = wavelet)):
                    coeffs = dtcwt(data = self.array, level = level, first_stage = first_stage, wavelet = wavelet)
                    reconstructed = idtcwt(coeffs = coeffs, first_stage = first_stage, wavelet = wavelet)
                    self.assertTrue(np.allclose(self.array, reconstructed))
    
    def test_axis(self):
        """ Test perfect reconstruction along all axes """
        for axis in range(0, self.array.ndim):
            coeffs = dtcwt(data = self.array, level = 2, axis = axis, first_stage = 'sym6', wavelet = 'qshift1')
            reconstructed = idtcwt(coeffs = coeffs, axis = axis, first_stage = 'sym6', wavelet = 'qshift1')
            self.assertTrue(np.allclose(self.array, reconstructed))
    
    def test_axis_limits(self):
        """ Test that an exception is raised for an invalid 'axis' parameter """
        with self.assertRaises(ValueError):
            coeffs = dtcwt(data = self.array, level = 1, axis = self.array.ndim, first_stage = 'sym6', wavelet = 'qshift1')
    
    def test_even_length_along_axis(self):
        """ Test that an exception is raised when array is not even along transform axis """
        with self.assertRaises(ValueError):
            dtcwt(data = np.zeros((17, 8)), level = 1, first_stage = 'sym6', wavelet = 'qshift1', axis = 0)

# Actual tests

class Test1D(TestDtcwt, unittest.TestCase):
    def setUp(self):
        self.array = np.random.random(size = (100,))

class Test2D(TestDtcwt, unittest.TestCase):
    def setUp(self):
        self.array = np.random.random(size = (50,50))

class Test3D(TestDtcwt, unittest.TestCase):
    def setUp(self):
        self.array = np.random.random(size = (10,10,10))
    
if __name__ == '__main__':
    unittest.main()