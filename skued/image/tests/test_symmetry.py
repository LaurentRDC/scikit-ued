# -*- coding: utf-8 -*-
import numpy as np
from .. import nfold
import unittest

np.random.seed(23)

class TestNFoldSymmetry(unittest.TestCase):

	def test_trivial(self):
		""" Test nfold_symmetry averaging on trivial array """
		im = np.zeros( (256, 256) )
		rot = nfold(im, center = (128, 128), mod = 3)
		self.assertTrue(np.allclose(rot, im))

	# TODO: test preserve_range = False has not effect

	def test_valid_mod(self):
		""" Test the the N-fold symmetry argument is valid """
		im = np.empty( (128, 128) )
		with self.assertRaises(ValueError):
			nfold(im, center = (64,64), mod = 1.7)
	
	def test_mask(self):
		""" Test that nfold_symmetry() works correctly with a mask """
		im = np.zeros((128, 128), dtype = np.int)
		mask = np.zeros_like(im, dtype = np.bool)

		im[0:20] = 1
		mask[0:20] = True
		
		rot = nfold(im, center = (64,64), mod = 2, mask = mask)
		self.assertTrue(np.allclose(rot, np.zeros_like(rot)))
	
	def test_no_side_effects(self):
		""" Test that nfold() does not modify the input image and mask """
		im = np.empty((128, 128), dtype = np.float)
		mask = np.zeros_like(im, dtype = np.bool)

		im.setflags(write = False)
		mask.setflags(write = False)

		rot = nfold(im, center = (67, 93),mod = 3, mask = mask)

if __name__ == '__main__':
	unittest.main()