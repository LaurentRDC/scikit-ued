# -*- coding: utf-8 -*-
import numpy as np
import unittest
from skimage import data

from .. import mnxc2

np.random.seed(23)
    
class TestMNXC2(unittest.TestCase):
	
	def test_side_effects(self):
		""" Test that mnxc2 does not modify the input in-place """
		im1 = np.random.random(size = (32,32,8))
		im2 = np.random.random(size = (32,32,8))

		m1 = np.random.choice([True, False], size = im1.shape)
		m2 = np.random.choice([True, False], size = im2.shape)

		# If arrays are written to, ValueError is raised
		for arr in (im1, im2, m1, m2):
			arr.setflags(write = False)
		
		xcorr = mnxc2(im1, im2, m1, m2)
	
	def test_range(self):
		im = np.random.random(size = (128,32,8))
		im2 = np.random.random(size = (128,32,8))

		m1 = np.random.choice([True, False], size = im.shape)
		m2 = np.random.choice([True, False], size = im2.shape)

		xcorr = mnxc2(im, im2, m1, m2).real
		self.assertGreaterEqual(xcorr.min(), -1)
		self.assertLessEqual(xcorr.max(), 1)
	
	def test_axes(self):
		""" Test that mnxc2 over axes is the same as a loop """
		im = np.random.random(size = (128,32,8))
		im2 = np.random.random(size = (128,32,8))

		m1 = np.random.choice([True, False], size = im.shape)
		m2 = np.random.choice([True, False], size = im2.shape)

		with_loop = np.empty_like(im)

		for index in range(im.shape[-1]):
			with_loop[:,:,index] = mnxc2(im[:,:,index], im2[:,:,index], m1[:,:,index], m2[:,:,index], axes = (0, 1), mode = 'same')
		
		over_axes = mnxc2(im, im2, m1, m2, axes = (0, 1), mode = 'same')

		self.assertTrue(np.allclose(with_loop, over_axes))


if __name__ == '__main__':
	unittest.main()