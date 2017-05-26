# -*- coding: utf-8 -*-
import numpy as np
from random import randint
import unittest
from skimage import data
from skimage.filters import gaussian
from skimage.transform import rotate

import matplotlib.pyplot as plt

from .. import shift_image, align
    
class TestAlign(unittest.TestCase):

	def test_trivial(self):
		""" Test alignment of identical images """
		aligned = tuple(align([data.camera() for _ in range(5)]))
		
		self.assertEqual(len(aligned), 5)
		self.assertSequenceEqual(data.camera().shape, aligned[0].shape)

	def test_misaligned_uniform_images(self):
		""" shift uniform images by entire pixels """
		original = np.ones((256, 256))
		misaligned = [shift_image(original, (1,-3))]
		aligned = align(misaligned, reference = original)

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

		aligned = align(misaligned, reference = original)

		# TODO: find a better figure-of-merit for alignment
		for im in aligned:
			# edge will be filled with zeros, we ignore
			diff = np.abs(original[5:-5, 5:-5] - im[5:-5, 5:-5])

			# Want less than 1% difference
			percent_diff = np.sum(diff) / (diff.size * (original.max() - original.min()))
			self.assertLess(percent_diff, 1)

if __name__ == '__main__':
	unittest.main()