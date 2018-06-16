# -*- coding: utf-8 -*-
from .. import powdersim
from ... import Crystal
from copy import deepcopy
import numpy as np
import unittest

class TestPowdersim(unittest.TestCase):

	def setUp(self):
		self.q = np.linspace(2, 10, 200)
		self.crystal = Crystal.from_database('C')

	def test_return_shape(self):
		""" Test that the return shape of powdersim() is as expected """
		pattern = powdersim(self.crystal, self.q)
		self.assertSequenceEqual(pattern.shape, self.q.shape)

if __name__ == '__main__':
	unittest.main()