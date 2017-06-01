# -*- coding: utf-8 -*-
import unittest
import numpy as np
from .. import repeated_array

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

if __name__ == '__main__':
	unittest.main()