# -*- coding: utf-8 -*-
from ..algorithms import (baseline_dt, baseline_dwt, _dwt_approx_rec, _dwt_approx_rec2, _dt_approx_rec)

import numpy as np
import unittest

np.random.seed(23)

class TestBaselineDWT(unittest.TestCase):
	def setUp(self):
		self.arr = np.random.random(size = (101,))

	def test_zero_level(self):
		""" baseline computed at the zero-th level should not affect the array """
		self.assertTrue(np.allclose(self.arr, baseline_dwt(self.arr, max_iter = 5, level = 0)))

	def test_approx_rec(self):
		""" Test that the underlying _dwt_approx_rec function is working properly """

		with self.subTest('1D shape'):
			arr = np.random.random(size = (102,))
			rec_arr = _dwt_approx_rec(arr, level = 2, wavelet = 'db6', mode = 'constant', axis = -1)
			self.assertSequenceEqual(rec_arr.shape, arr.shape)

		with self.subTest('1D along axis'):
			arr2 = np.random.random(size = (101, 104))
			rec_arr2 = _dwt_approx_rec(arr2, level = 2, wavelet = 'sym4', mode = 'constant', axis = 1)
			self.assertSequenceEqual(rec_arr2.shape, arr2.shape)
		
		with self.subTest('2D shape'):
			arr2 = np.random.random(size = (102,104))
			rec_arr2 = _dwt_approx_rec2(arr2, level = 2, wavelet = 'sym6', mode = 'constant', axis = (-2, -1))
			self.assertSequenceEqual(rec_arr2.shape, arr2.shape)

		with self.subTest('2D along axes'):
			arr2 = np.random.random(size = (102,94,25))
			rec_arr2 = _dwt_approx_rec2(arr2, level = 2, wavelet = 'sym6', mode = 'constant', axis = (0, 1))
			self.assertSequenceEqual(rec_arr2.shape, arr2.shape)

	def test_trivial_case_1d(self):
		""" The baseline for a 1d array of zeros should be zeros """
		arr = np.zeros_like(self.arr)
		self.assertTrue(np.allclose(arr, baseline_dwt(arr, max_iter = 10, level = None)))

	def test_trivial_case_2d(self):
		""" The baseline for a 2d array of zeros should be zeros """
		arr = np.zeros( (21, 31, 5))
		self.assertTrue(np.allclose(arr, baseline_dwt(arr, max_iter = 10, level = None, axis = (0, 1))))

	def test_1d_along_axis(self):
		""" Test that iterating over array rows and baseline_dwt along axis are equivalent """
		block = np.random.random(size = (21, 51))

		# Iterate over rows
		baseline_iterated = np.empty_like(block)
		for index, row in enumerate(block):
			baseline_iterated[index, :] = baseline_dwt(row, max_iter = 50)

		# along axis
		baseline_axis = baseline_dwt(block, max_iter = 50, axis = 1)

		self.assertTrue(np.allclose(baseline_axis, baseline_iterated))

class TestbaselineDT(unittest.TestCase):
	def setUp(self):
		self.arr = np.random.random(size = (101,))

	def test_zero_level(self):
		""" baseline computed at the zero-th level should not affect the array """
		self.assertTrue(np.allclose(self.arr, baseline_dt(self.arr, max_iter = 5, level = 0)))

	def test_approx_rec(self):
		""" Test that the underlying _dwt_approx_rec function is working properly """

		with self.subTest('1D shape'):
			arr = np.random.random(size = (102,))
			rec_arr = _dt_approx_rec(arr, level = 2, first_stage = 'db1', wavelet = 'qshift3', mode = 'smooth', axis = -1)
			self.assertSequenceEqual(rec_arr.shape, arr.shape)

		with self.subTest('2D along axis'):
			arr2 = np.random.random(size = (21, 52))
			rec_arr2 = _dt_approx_rec(arr2, level = 2, first_stage = 'db1', wavelet = 'qshift3', mode = 'smooth', axis = 1)
			self.assertSequenceEqual(rec_arr2.shape, arr2.shape)

	def test_trivial_case(self):
		""" The baseline for an array of zeros should be zeros """
		arr = np.zeros_like(self.arr)
		self.assertTrue(np.allclose(arr, baseline_dt(arr, max_iter = 10, level = None)))
	
	def test_positive_baseline(self):
		""" Test that the baseline is never negative """
		arr = 10*np.random.random(size = (128,))
		baseline = baseline_dt(arr, max_iter = 10)

		self.assertTrue(np.all(np.greater_equal(baseline, 0)))

	def test_baseline_limit(self):
		""" Test that the baseline is never more than the original signal """
		arr = 10*np.random.random(size = (128,))
		baseline = baseline_dt(arr, max_iter = 10)
		
		self.assertTrue(np.all(np.greater_equal(arr, baseline)))
	
	def test_final_shape(self):
		""" Test that baseline_dt returns an array of the same shape as input """
		
		with self.subTest('1D baseline of odd length'):
			arr = np.random.random(size = (101,))
			b = baseline_dt(arr, max_iter = 10, level = None)
			self.assertSequenceEqual(arr.shape, b.shape)
		
		with self.subTest('1D baseline of even length'):
			arr = np.random.random(size = (100,))
			b = baseline_dt(arr, max_iter = 10, level = None)
			self.assertSequenceEqual(arr.shape, b.shape)
		
		with self.subTest('baseline along axis of odd length'):
			arr = np.random.random(size = (101, 113))
			b = baseline_dt(arr, max_iter = 10, level = None, axis = 0)
			self.assertSequenceEqual(arr.shape, b.shape)

		with self.subTest('baseline along axis of even length'):
			arr = np.random.random(size = (102, 104))
			b = baseline_dt(arr, max_iter = 10, level = None, axis = 1)
			self.assertSequenceEqual(arr.shape, b.shape)

	def test_2d_along_axis(self):
		""" Test that iterating over array rows and baseline_dt along axis are equivalent """
		block = np.random.random(size = (21, 51))

		# Iterate over rows
		baseline_iterated = np.empty_like(block)
		for index, row in enumerate(block):
			baseline_iterated[index, :] = baseline_dt(row, max_iter = 50)

		# along axis
		baseline_axis = baseline_dt(block, max_iter = 50, axis = 1)

		self.assertTrue(np.allclose(baseline_axis, baseline_iterated))
	
	def test_background_regions(self):
		""" Test that background_regions is used correctly along certain axes """
		block = np.random.random(size = (21, 51))

		baseline_params = {'max_iter': 50,
						   'background_regions': [(..., range(0, 3))]}

		# Iterate over rows
		baseline_iterated = np.empty_like(block)
		for index, row in enumerate(block):
			baseline_iterated[index, :] = baseline_dt(row, **baseline_params)

		# along axis
		baseline_axis = baseline_dt(block, axis = 1, **baseline_params)

		self.assertTrue(np.allclose(baseline_axis, baseline_iterated))


if __name__ == '__main__':
    unittest.main()



