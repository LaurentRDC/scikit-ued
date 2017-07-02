# -*- coding: utf-8 -*-
import numpy as np
from random import randint
import unittest
from skimage import data
from skimage.filters import gaussian
from skimage.transform import rotate

from .. import shift_image, align, diff_register
from .test_powder import circle_image

class TestDiffRegister(unittest.TestCase):
	""" Test the diff_register function """

	def test_on_skimage_data(self):
		""" Test registering translation from scikit-image's data module """
		with self.subTest('skimage.data.camera()'):
			im = data.camera()
			shifted = shift_image(im, (4, -5))

			shift = diff_register(shifted, reference = im, search_space = 7)
			self.assertSequenceEqual(tuple(shift), (4, -5))
		
		with self.subTest('skimage.data.coins()'):
			im = data.coins()
			shifted = shift_image(im, (0, 3))

			shift = diff_register(shifted, reference = im, search_space = 3)
			self.assertSequenceEqual(tuple(shift), (0, 3))
	
	def test_on_simulated_powder(self):
		""" Test on simulated perfect powder rings """
		center = (64, 64)
		im = circle_image(shape = (128, 128), center = center, 
						  radii = [16, 32], intensities = [2,1])
		im[:] = gaussian(im, 2)

		shifted = shift_image(im, (-2, 6))
		shift = diff_register(shifted, reference = im, search_space = 6)
		self.assertSequenceEqual(tuple(shift), (-2, 6))
    
	def test_on_masked_skimage(self):
		""" Test diff_register on masked image data from scikit-image """
		im = np.asfarray(data.camera())
		shifted = shift_image(im, (5, -2))

		mask = np.zeros_like(im, dtype = np.bool)
		im[256, 256] = 1e10
		shifted[256, 256] = 1e10
		mask[256, 256] = True

		shift = diff_register(shifted, reference = im, mask = mask, search_space = 6)
		self.assertSequenceEqual(tuple(shift), (5, -2))

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