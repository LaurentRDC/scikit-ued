# -*- coding: utf-8 -*-
import numpy as np
import unittest
from skimage import data

from .. import mnxc2, register_translation, shift_image

np.random.seed(23)
    
class TestMNXC2(unittest.TestCase):

	def test_trivial(self):
		""" Test that the autocorrelation using mnxc2 has a single
		peak of 1 """
		im = np.random.random(size = (32,32))
		im2 = np.random.random(size = (32,32))

		xcorr = mnxc2(im, im2)
		self.assertGreaterEqual(xcorr.min(), -1)
		self.assertLessEqual(xcorr.max(), 1)
	
	def test_side_effects(self):
		""" Test that mnxc2 does not modify the input in-place """
		im1 = np.random.random(size = (32,32))
		im2 = np.random.random(size = (32,32))

		m1 = np.random.choice([True, False], size = im1.shape)
		m2 = np.random.choice([True, False], size = im2.shape)

		# If arrays are written to, ValueError is raised
		for arr in (im1, im2, m1, m2):
			arr.setflags(write = False)
		
		xcorr = mnxc2(im1, im2, m1, m2)
	
	def test_autocorr(self):
		""" Test mnxc2 on two identical images for a peak in the center """
		im = np.random.random(size = (128, 128))

		with self.subTest('No masks'):
			xcorr = mnxc2(im, im, mode = 'same')
			center = np.unravel_index(np.argmax(xcorr), xcorr.shape)

			self.assertAlmostEqual(xcorr[center], 1)
			# Due to the way arrays are centered in mnxc2, center is at [63,63]
			self.assertTrue(np.allclose(center,np.array(im.shape)/2 - 1, atol = 1))
	
	def test_range(self):
		im = np.random.random(size = (128,32))
		im2 = np.random.random(size = (128,32))

		m1 = np.random.choice([True, False], size = im.shape)
		m2 = np.random.choice([True, False], size = im2.shape)

		xcorr = mnxc2(im, im2, m1, m2).real
		self.assertGreaterEqual(xcorr.min(), -1)
		self.assertLessEqual(xcorr.max(), 1)

class TestRegisterTranslation(unittest.TestCase):

	def test_trivial_skimage_data(self):
		""" Test that the translation between two identical images is (0,0), even
		with added random noise and random masks """

		im = np.asfarray(data.camera())

		with self.subTest('No noise'):
			shift = register_translation(im, im)

			self.assertTrue(np.allclose(shift, (0,0), atol = 1))
		
		with self.subTest('With 5% noise'):
			noise1 = 0.05 * im.max() * np.random.random(size = im.shape)
			noise2 = 0.05 * im.max() * np.random.random(size = im.shape)

			shift = register_translation(im + noise1, im + noise2)
			self.assertTrue(np.allclose(shift, (0,0), atol = 1))
		
		with self.subTest('With random masks'):
			m1 = np.random.choice([True, False], size = im.shape)
			m2 = np.random.choice([True, False], size = im.shape)

			shift = register_translation(im, im, m1, m2)
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
			shift = register_translation(im, im2, edge_mask)

			self.assertTrue(np.allclose(shift, -random_shift, atol = 1))
		
		with self.subTest('With 5% noise'):
			noise1 = 0.05 * im.max() * np.random.random(size = im.shape)
			noise2 = 0.05 * im.max() * np.random.random(size = im.shape)

			shift = register_translation(im + noise1, im2 + noise2, edge_mask)
			self.assertTrue(np.allclose(shift, -random_shift, atol = 1))
		
		with self.subTest('With random masks'):
			m1 = np.random.choice([True, False], size = im.shape)
			m2 = np.random.choice([True, False], size = im.shape)

			shift = register_translation(im, im2, m1, m2)
			self.assertTrue(np.allclose(shift, -random_shift, atol = 1))


if __name__ == '__main__':
	unittest.main()