from ..algorithms import (baseline_dt, baseline_dwt)

import numpy as np
import unittest

np.random.seed(23)

class TestBaselineDWT(unittest.TestCase):
	def setUp(self):
		self.arr = np.random.random(size = (101,))

	def test_zero_level(self):
		""" baseline computed at the zero-th level should not affect the array """
		self.assertTrue(np.allclose(self.arr, baseline_dwt(self.arr, max_iter = 5, level = 0)))

	def test_trivial_case(self):
		""" The baseline for an array of zeros should be zeros """
		arr = np.zeros_like(self.arr)
		self.assertTrue(np.allclose(arr, baseline_dwt(arr, max_iter = 10, level = 'max')))

class TestbaselineDT(unittest.TestCase):
	def setUp(self):
		self.arr = np.random.random(size = (101,))

	def test_zero_level(self):
		""" baseline computed at the zero-th level should not affect the array """
		self.assertTrue(np.allclose(self.arr, baseline_dt(self.arr, max_iter = 5, level = 0)))

	def test_trivial_case(self):
		""" The baseline for an array of zeros should be zeros """
		arr = np.zeros_like(self.arr)
		self.assertTrue(np.allclose(arr, baseline_dt(arr, max_iter = 10, level = 'max')))

if __name__ == '__main__':
    unittest.main()



