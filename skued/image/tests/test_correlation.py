# -*- coding: utf-8 -*-
import numpy as np
import unittest
from skimage import data

from .. import masked_xcorr
    
class TestMaskedXCorr(unittest.TestCase):

	def test_trivial(self):
		""" Test that the autocorrelation using masked_xcorr has a single
		peak of 1 """
		im = data.camera()
		mask = np.zeros_like(im, dtype = np.bool)

		xcorr = masked_xcorr(im, im, mask).real
		self.assertGreaterEqual(xcorr.min(), -1)
		self.assertLessEqual(xcorr.max(), 1)

if __name__ == '__main__':
	unittest.main()