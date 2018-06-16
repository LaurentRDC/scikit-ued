# -*- coding: utf-8 -*-
import unittest
import numpy as np
from .. import (repeated_array, mirror, cart2polar, polar2cart, plane_mesh,
                spherical2cart, cart2spherical, complex_array)

np.random.seed(23)

class TestRepeatedArray(unittest.TestCase):

	def setUp(self):
		self.arr = np.random.random(size = (4,5))

	def test_trivial(self):
		""" Test repeated_array of 0 copies """
		composite = repeated_array(self.arr, num = 0, axes = 0)
		self.assertTrue(np.allclose(composite, self.arr))
	
	def test_single_axis(self):
		""" Test repeating the array over a single axis """
		composite = repeated_array(self.arr, num = 3, axes = 1)
		expected_new_shape = (self.arr.shape[0], self.arr.shape[1]*3)
		self.assertEqual(composite.shape, expected_new_shape)
	
	def test_multiple_axes(self):
		""" Test repeating an array over multiple axes in all possible orders """
		with self.subTest('axes = (0, 1)'):
			composite = repeated_array(self.arr, num = (3, 2), axes = (0, 1))
			expected_new_shape = (self.arr.shape[0]*3, self.arr.shape[1]*2)
			self.assertEqual(composite.shape, expected_new_shape)

		with self.subTest('axes = (1, 0)'):
			composite = repeated_array(self.arr, num = (2, 3), axes = (1, 0))
			expected_new_shape = (self.arr.shape[0]*3, self.arr.shape[1]*2)
			self.assertEqual(composite.shape, expected_new_shape)

class TestComplexArray(unittest.TestCase):

	def test_floats(self):
		""" Test that two floating arrays are cast correctly """
		real, imag = np.empty((3,4), dtype = np.float), np.empty((3,4), dtype = np.float)
		self.assertEqual(complex_array(real, imag).dtype, np.complex)
	
	def test_non_floats(self):
		""" Test that two integer arrays are cast correctly """
		real, imag = np.empty((3,4), dtype = np.int16), np.empty((3,4), dtype = np.int8)
		self.assertEqual(complex_array(real, imag).dtype, np.complex)
	
	def test_results(self):
		""" test that ``complex_array`` returns appropriate results """
		arr1 = np.random.random((4, 5))
		arr2 = np.random.random((4, 5))

		from_complex_array = complex_array(arr1, arr2)
		by_hand = arr1 + 1j*arr2
		self.assertEqual(from_complex_array.dtype, by_hand.dtype)
		self.assertTrue(np.allclose(from_complex_array, by_hand))

class TestMirror(unittest.TestCase):

	def test_1D(self):
		""" Test mirror() on a 1D array """
		arr = np.zeros( (16,), dtype = np.float)
		arr[15] = 1
		self.assertTrue(np.allclose(arr[::-1], mirror(arr)))

	def test_2D_all_axes(self):
		""" Test mirror() on a 2D array for all axes """
		arr = np.zeros( (16,16), dtype = np.float)
		arr[15, 3] = 1
		self.assertTrue(np.allclose(arr[::-1, ::-1], mirror(arr)))
	
	def test_2D_one_axis(self):
		""" Test mirror() on a 2D array for one axis """
		arr = np.zeros( (16,16), dtype = np.float)
		arr[15, 3] = 1
		self.assertTrue(np.allclose(arr[:, ::-1], mirror(arr, axes = 1)))
		self.assertTrue(np.allclose(arr[::-1, :], mirror(arr, axes = 0)))

class TestCart2Polar(unittest.TestCase):

    def test_back_and_forth(self):
        """ Test that cart2polar and polar2cart are reciprocal """
        x = np.random.random(size = (16, 8))
        y = np.random.random(size = (16, 8))

        r, t = cart2polar(x, y)

        xp, yp = polar2cart(r,t)
        self.assertTrue(np.allclose(x, xp))
        self.assertTrue(np.allclose(y, yp))

class TestSpherical2Cart(unittest.TestCase):

    def test_back_and_forth(self):
        """ Test that cart2polar and polar2cart are reciprocal """
        x = np.random.random(size = (16, 8))
        y = np.random.random(size = (16, 8))
        z = np.random.random(size = (16, 8))

        r, p, t = cart2spherical(x, y, z)

        xp, yp, zp = spherical2cart(r,p, t)
        self.assertTrue(np.allclose(x, xp))
        self.assertTrue(np.allclose(y, yp))
        self.assertTrue(np.allclose(z, zp))
        
class TestPlaneMesh(unittest.TestCase):

    def test_shape(self):
        """ Test that shape is as expected """
        extent1 = np.linspace(0, 10, num = 64)
        extent2 = np.linspace(0, 10, num = 128)
        v1, v2, _ = np.eye(3)
        
        for arr in plane_mesh(v1, v2, extent1, extent2):
            self.assertSequenceEqual(arr.shape, (64, 128))
    
    def test_origin(self):
        """ Test that plane_mesh is generated from origin """
        extent1 = np.linspace(0, 10, num = 64)
        extent2 = np.linspace(0, 10, num = 128)
        v1, v2, _ = np.eye(3)

        for arr in plane_mesh(v1, v2, extent1, extent2, origin = (-4, -4, -4)):
            self.assertEqual(arr.min(), -4)
        
		
if __name__ == '__main__':
	unittest.main()