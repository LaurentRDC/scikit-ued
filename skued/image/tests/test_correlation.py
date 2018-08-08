# -*- coding: utf-8 -*-
import unittest
from pathlib import Path

import numpy as np
from scipy.signal import correlate
from skimage.io import imread

from .. import mnxc2, xcorr

np.random.seed(23)

def _masked_register_translation(im1, im2, m1, m2, overlap_ratio):
	xcorr = mnxc2(im1, im2, m1, m2, mode = 'full', overlap_ratio = overlap_ratio)
	# Generalize to the average of multiple maxima
	maxima = np.transpose(np.nonzero(xcorr == xcorr.max()))
	center = np.mean(maxima, axis = 0)
	shift = center - np.array(im2.shape) + 1
	return -shift

class TestXcorr(unittest.TestCase):

	def test_final_shape(self):
		""" Test that the full cross correlation has shape s1 + s2  - 1 """
		im1 = np.random.random(size = (32,32,5))
		im2 = np.random.random(size = (32,32,5))
		ret = xcorr(im1, im2, mode = 'full', axes = (0,1,2))
		self.assertTupleEqual(ret.shape, (63,63,9))
	
	def test_against_scipy_float(self):
		""" Test against scipy.signal.correlate with float inputs """
		im1 = np.random.random(size = (32,32,5))
		im2 = np.random.random(size = (32,32,5))

		with self.subTest("mode = 'full'"):
			from_skued = xcorr(im1, im2, mode = 'full', axes = (0,1,2))
			from_scipy = correlate(im1, im2, mode = 'full')

			self.assertTrue(np.allclose(from_scipy, from_skued))

		with self.subTest("mode = 'same'"):
			from_skued = xcorr(im1, im2, mode = 'same', axes = (0,1,2))
			from_scipy = correlate(im1, im2, mode = 'same')

			self.assertTrue(np.allclose(from_scipy, from_skued))

	def test_against_scipy_complex(self):
		""" Test against scipy.signal.correlate with complex inputs """
		im1 = np.random.random(size = (32,32,5)) + 1j*np.random.random(size = (32,32,5))
		im2 = np.random.random(size = (32,32,5)) + 1j*np.random.random(size = (32,32,5))

		with self.subTest("mode = 'full'"):
			from_skued = xcorr(im1, im2, mode = 'full', axes = (0,1,2))
			from_scipy = correlate(im1, im2, mode = 'full')

			self.assertTrue(np.allclose(from_scipy, from_skued))

		with self.subTest("mode = 'same'"):
			from_skued = xcorr(im1, im2, mode = 'same', axes = (0,1,2))
			from_scipy = correlate(im1, im2, mode = 'same')
			
			self.assertTrue(np.allclose(from_scipy, from_skued))
    
class TestMNXC2(unittest.TestCase):
	
	def test_side_effects(self):
		""" Test that mnxc2 does not modify the input in-place """
		im1 = np.random.random(size = (32,32,5))
		im2 = np.random.random(size = (32,32,5))

		m1 = np.random.choice([True, False], size = im1.shape)
		m2 = np.random.choice([True, False], size = im2.shape)

		# If arrays are written to, ValueError is raised
		for arr in (im1, im2, m1, m2):
			arr.setflags(write = False)
		
		xcorr = mnxc2(im1, im2, m1, m2)
	
	def test_range(self):
		im = np.random.random(size = (128,32,5))
		im2 = np.random.random(size = (128,32,5))

		m1 = np.random.choice([True, False], size = im.shape)
		m2 = np.random.choice([True, False], size = im2.shape)

		xcorr = mnxc2(im, im2, m1, m2).real
		self.assertGreaterEqual(xcorr.min(), -1)
		self.assertLessEqual(xcorr.max(), 1)
	
	def test_axes(self):
		""" Test that mnxc2 over axes is the same as a loop """
		im = np.random.random(size = (128,32,5))
		im2 = np.random.random(size = (128,32,5))

		m1 = np.random.choice([True, False], size = im.shape)
		m2 = np.random.choice([True, False], size = im2.shape)

		with_loop = np.empty_like(im)

		for index in range(im.shape[-1]):
			with_loop[:,:,index] = mnxc2(im[:,:,index], im2[:,:,index], m1[:,:,index], m2[:,:,index], axes = (0, 1), mode = 'same')
		
		over_axes = mnxc2(im, im2, m1, m2, axes = (0, 1), mode = 'same')

		self.assertTrue(np.allclose(with_loop, over_axes))

	def test_final_shape(self):
		""" Test that the MNXC2 has shape s1 + s2  - 1 """
		im1 = np.random.random(size = (32,32,5))
		im2 = np.random.random(size = (32,32,5))
		ret = mnxc2(im1, im2)
		self.assertTupleEqual(ret.shape, (63,63,5))

	def test_padfield_data(self):
		""" Test translation registration for data included in Padfield 2010 """
		# Test translated from MATLABimplementation `MaskedFFTRegistrationTest` file. You can find the source code here: 
		# http://www.dirkpadfield.com/Home/MaskedFFTRegistrationCode.zip
		IMAGES_DIR = Path(__file__).parent / 'images'

		shifts = [(75, 75), (-130, 130), (130, 130)]
		for xi, yi in shifts:
			with self.subTest('X = {:d}, Y = {:d}'.format(xi, yi)):
				fixed_image = imread(IMAGES_DIR / 'OriginalX{:d}Y{:d}.png'.format(xi, yi))
				moving_image = imread(IMAGES_DIR/ 'TransformedX{:d}Y{:d}.png'.format(xi, yi))

				# Our definition for masks is inverted from Padfields
				# Invalid pixels are 1
				fixed_mask = (fixed_image == 0)
				moving_mask = (moving_image == 0)

				# Note that shifts in x and y and shifts in cols and rows
				shift_y, shift_x = _masked_register_translation(fixed_image, moving_image, fixed_mask, moving_mask, overlap_ratio = 1/10)
				# NOTE: by looking at the test code from Padfield's MaskedFFTRegistrationCode repository,
				#		the shifts were not xi and yi, but xi and -yi
				self.assertTupleEqual((xi, -yi), (shift_x, shift_y))

if __name__ == '__main__':
	unittest.main()
