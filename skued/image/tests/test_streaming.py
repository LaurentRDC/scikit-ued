# -*- coding: utf-8 -*-
import unittest
from itertools import repeat
from random import randint, random

import numpy as np
from scipy.stats import sem as scipy_sem
from skimage import data

from .. import ialign, iaverage, isem, shift_image, last


class TestIAlign(unittest.TestCase):

	def test_trivial(self):
		""" Test alignment of identical images """
		aligned = tuple(ialign([data.camera() for _ in range(5)]))
		
		self.assertEqual(len(aligned), 5)
		self.assertSequenceEqual(data.camera().shape, aligned[0].shape)

	def test_misaligned_uniform_images(self):
		""" shift uniform images by entire pixels """
		original = np.ones((256, 256))
		misaligned = [shift_image(original, (1,-3))]
		aligned = ialign(misaligned, reference = original)

		for im in aligned:
			# edge will be filled with zeros
			self.assertTrue(np.allclose(original[5:-5, 5:-5], 
										im[5:-5, 5:-5]))

	def test_misaligned_canned_images(self):
		""" shift images from skimage.data by entire pixels.
	   We don't expect perfect alignment."""
		original = data.camera()
		misaligned = [shift_image(original, (randint(-4, 4), randint(-4, 4))) 
					  for _ in range(5)]

		aligned = ialign(misaligned, reference = original)

		# TODO: find a better figure-of-merit for alignment
		for im in aligned:
			# edge will be filled with zeros, we ignore
			diff = np.abs(original[5:-5, 5:-5] - im[5:-5, 5:-5])

			# Want less than 1% difference
			percent_diff = np.sum(diff) / (diff.size * (original.max() - original.min()))
			self.assertLess(percent_diff, 1)

class TestIAverage(unittest.TestCase):

	def test_trivial(self):
		""" Test iaverage on stream of zeroes """
		stream = repeat(np.zeros( (64,64), dtype = np.float ), times = 5)
		for av in iaverage(stream):
			self.assertTrue(np.allclose(av, np.zeros_like(av)))
	
	def test_mean(self):
		""" Test results against of unweighted average against numpy.mean """
		stream = [np.random.random(size = (64,64)) for _ in range(5)]

		from_iaverage = last(iaverage(stream))
		from_numpy = np.mean(np.dstack(stream), axis = 2)

		self.assertTrue(np.allclose(from_iaverage, from_numpy))
	
	def test_average(self):
		""" Test results of weighted average against numpy.average """
		stream = [np.random.random(size = (2,2)) for _ in range(5)]

		with self.subTest('float weights'):
			weights = [random() for _ in stream]
			from_iaverage = last(iaverage(stream, weights = weights))
			from_numpy = np.average(np.dstack(stream), axis = 2, weights = np.array(weights))

			self.assertTrue(np.allclose(from_iaverage, from_numpy))
		
		with self.subTest('array weights'):
			weights = [np.random.random(size = stream[0].shape) for _ in stream]
			from_iaverage = last(iaverage(stream, weights = weights))
			from_numpy= np.average(np.dstack(stream), axis = 2, weights = np.dstack(weights))

			self.assertTrue(np.allclose(from_iaverage, from_numpy))

class TestISem(unittest.TestCase):

	def test_first(self):
		""" Test that the first yielded value of isem is an array fo zeros """
		stream = repeat(np.random.random( size = (64,64)), times = 5)
		first = next(isem(stream))

		self.assertTrue(np.allclose(first, np.zeros_like(first)))
	
	def test_against_scipy_sem(self):
		""" Test that the results of isem are in agreement with scipy.stats.sem """
		stream = [np.random.random(size = (64,64)) for _ in range(5)]

		from_isem = last(isem(stream))
		from_scipy = scipy_sem(np.dstack(stream), axis = 2, ddof = 1)

		self.assertTrue(np.allclose(from_isem, from_scipy))

if __name__ == '__main__':
	unittest.main()
