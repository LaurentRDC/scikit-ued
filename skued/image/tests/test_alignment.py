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

	def test_trivial_skimage_data(self):
		""" Test that the translation between two identical images is (0,0), even
		with added random noise and random masks """

		im = np.asfarray(data.camera())

		with self.subTest('No noise'):
			shift = diff_register(im, im)

			self.assertTrue(np.allclose(shift, (0,0), atol = 1))
		
		with self.subTest('With 5% noise'):
			noise1 = 0.05 * im.max() * np.random.random(size = im.shape)
			noise2 = 0.05 * im.max() * np.random.random(size = im.shape)

			shift = diff_register(im + noise1, im + noise2)
			self.assertTrue(np.allclose(shift, (0,0), atol = 1))
		
		with self.subTest('With random masks'):
			m1 = np.random.choice([True, False], size = im.shape)

			shift = diff_register(im, im, m1)
			self.assertTrue(np.allclose(shift, (0,0), atol = 1))
	
	def test_shifted_skimage_data(self):
		""" Test that translation is registered for data from scikit-image """
		random_shift = np.random.randint(low = 0, high = 10, size = (2,))

		im = np.asfarray(data.camera())
		im2 = np.asfarray(data.camera())
		im2[:] = shift_image(im2, shift = random_shift, fill_value = 0)

		# Masking the edges due to shifting
		edge_mask = np.ones_like(im, dtype = np.bool)
		edge_mask[6:-6, 6:-6] = False

		with self.subTest('No noise'):
			shift = diff_register(im, im2, edge_mask)

			self.assertTrue(np.allclose(shift, -random_shift, atol = 1))
		
		with self.subTest('With 5% noise'):
			noise1 = 0.05 * im.max() * np.random.random(size = im.shape)
			noise2 = 0.05 * im.max() * np.random.random(size = im.shape)

			shift = diff_register(im + noise1, im2 + noise2, edge_mask)
			self.assertTrue(np.allclose(shift, -random_shift, atol = 1))
		
		with self.subTest('With random mask'):
			m1 = np.random.choice([True, False], size = im.shape)

			shift = diff_register(im, im2, m1)
			self.assertTrue(np.allclose(shift, -random_shift, atol = 1))
	
	def test_side_effects(self):
		""" Test that arrays registered by diff_register are not modified """
		im1 = np.random.random(size = (32,32))
		im2 = np.random.random(size = (32,32))
		mask = np.random.choice([True, False], size = im1.shape)

		# If arrays are written to, ValueError is raised
		for arr in (im1, im2, mask):
			arr.setflags(write = False)
		
		shift = diff_register(im1, im2, mask)

class TestAlign(unittest.TestCase):

	def test_trivial(self):
		""" Test alignment of identical images """
		aligned = align(data.camera(), reference = data.camera())
		self.assertTrue(np.allclose(aligned, data.camera()))

	def test_misaligned_uniform_images(self):
		""" shift uniform images by entire pixels """
		original = np.ones((256, 256))
		misaligned = shift_image(original, (1,-3))
		aligned = align(misaligned, reference = original)
		
		self.assertTrue(np.allclose(original[5:-5, 5:-5], 
									aligned[5:-5, 5:-5]))

	def test_misaligned_canned_images(self):
		""" shift images from skimage.data by entire pixels.
	   We don't expect perfect alignment."""
		original = data.camera()
		misaligned = shift_image(original, (randint(-4, 4), randint(-4, 4))) 

		aligned = align(misaligned, reference = original)

		# edge will be filled with zeros, we ignore
		diff = np.abs(original[5:-5, 5:-5] - aligned[5:-5, 5:-5])

		# Want less than 1% difference
		percent_diff = np.sum(diff) / (diff.size * (original.max() - original.min()))
		self.assertLess(percent_diff, 1)

if __name__ == '__main__':
	unittest.main()