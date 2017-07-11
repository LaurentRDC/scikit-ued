# -*- coding: utf-8 -*-
import numpy as np
from .. import nfold_symmetry
import unittest

np.random.seed(23)

class TestNFoldSymmetry(unittest.TestCase):

	def test_trivial(self):
		""" Test nfold_symmetry averaging on trivial array """
		im = np.zeros( (256, 256) )
		rot = nfold_symmetry(im, center = (128, 128), mod = 3)
		self.assertTrue(np.allclose(rot, im))

	# TODO: test preserve_range = False has not effect

	def test_valid_mod(self):
		""" Test the the N-fold symmetry argument is valid """
		im = np.empty( (128, 128) )
		with self.assertRaises(ValueError):
			nfold_symmetry(im, center = (64,64), mod = 1.7)

if __name__ == '__main__':
	unittest.main()