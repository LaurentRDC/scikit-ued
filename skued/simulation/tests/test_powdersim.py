# -*- coding: utf-8 -*-
from .. import powdersim
from ... import Crystal
from copy import deepcopy
import numpy as np
import unittest

class TestPowdersim(unittest.TestCase):

	def setUp(self):
		self.scattering_length = np.linspace(0.11, 0.8, 200)
		self.crystal = Crystal.from_database('C')

	def test_return_shape(self):
		""" Test that the return shape of powdersim() is as expected """
		pattern = powdersim(self.crystal, self.scattering_length)
		self.assertSequenceEqual(pattern.shape, self.scattering_length.shape)

if __name__ == '__main__':
	unittest.main()