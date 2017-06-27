# -*- coding: utf-8 -*-
import numpy as np
from random import randint
import unittest
from skimage import data
from skimage.filters import gaussian
from skimage.transform import rotate

from .. import shift_image, align, register_translation

np.random.seed(23)
    
class TestAlign(unittest.TestCase):

	def test_trivial(self):
		""" Test alignment of identical images """
		aligned = tuple(align([data.camera() for _ in range(5)]))
		
		self.assertEqual(len(aligned), 5)
		self.assertSequenceEqual(data.camera().shape, aligned[0].shape)

	def test_misaligned_uniform_images_no_mask(self):
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


class TestRegisterTranslation(unittest.TestCase):
	""" Test skued.image_analysis.register_translation """

	def test_trivial(self):
		""" Test that the translation between two identical images is zero """
		# If both masks were None, register_translation would fall back to skimage.feature.register_translation
		# Thus, we use a trivial mask

		im = data.camera()
		mask = np.ones_like(im, dtype = np.bool)

		shifts = register_translation(im, im, mask, mask)
		self.assertSequenceEqual((0,0), tuple(shifts))
	
	def test_shift_trivial_mask(self):
		""" Test that register_translation registers translations when using trivial masks
		(without falling back to scikit-image's implementation) """

		im = data.camera()
		im2 = shift_image(im, shift = (-3,1))

		mask = np.ones_like(im, dtype = np.bool)

		shifts = register_translation(im, im2, mask, mask)
		self.assertSequenceEqual((3, -1), tuple(shifts))

	def test_shift_random_mask(self):
		""" Test that register_translation registers translations when using random masks """

		im = data.camera()
		im2 = shift_image(im, shift = (-3,1))

		with self.subTest('Single mask'):

			mask = np.ones_like(im, dtype = np.bool)
			mask[400:512, 256:512] = False

			shifts = register_translation(im, im2, mask, mask)
			self.assertSequenceEqual((3, -1), tuple(shifts))
		
		with self.subTest('Two masks'):

			mask1 = np.ones_like(im, dtype = np.bool)
			mask2 = np.array(mask1, copy = True)

			mask1[400:512, 256:512] = False
			mask2[0:256, 0:128] = False

			shifts = register_translation(im, im2, mask1, mask2)
			self.assertSequenceEqual((3, -1), tuple(shifts))

	def test_shift_bright_pixels(self):
		""" Test that register_translation registers translations correctly in the presence
		of extremely bright fixed pixels that are masked. """

		im = data.camera()
		im2 = shift_image(im, shift = (-3,1))

		im[258:260, 258:260] = 1e25
		im2[258:260, 258:260] = 1e25
		
		mask = np.ones_like(im)
		mask[256:262, 256:262] = False

		shifts = register_translation(im, im2, mask, mask)
		self.assertSequenceEqual((3, -1), tuple(shifts))

	def test_shift_no_mask(self):
		""" Test that register_translation registers translations even with no masks,
		in which case it falls back to skimage.feature.register_translation """

		im = data.camera()
		im2 = shift_image(im, shift = (-3,1))

		mask = np.ones_like(im, dtype = np.bool)

		shifts = register_translation(im, im2, None, None)
		self.assertSequenceEqual((3, -1), tuple(shifts))

if __name__ == '__main__':
	unittest.main()