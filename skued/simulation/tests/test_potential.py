# -*- coding: utf-8 -*-
from .. import electrostatic, pelectrostatic
from ...structure import graphite
from copy import deepcopy
import numpy as np
import unittest

class TestElectrostatic(unittest.TestCase):

	def setUp(self):
		self.crystal = deepcopy(graphite)

	def test_return_shape(self):
		""" Test that the return shape of pelectrostatic is the same as input arrays """
		xx, yy, zz = np.meshgrid(np.linspace(-10, 10, 16), np.linspace(-10, 10, 16), np.linspace(-10, 10, 16))
		potential = electrostatic(graphite, xx, yy, zz)

		self.assertSequenceEqual(xx.shape, potential.shape)
	
	def test_side_effects(self):
		""" Test that mesh arrays are not written to in pelectrostatic """
		xx, yy, zz = np.meshgrid(np.linspace(-10, 10, 16), np.linspace(-10, 10, 16), np.linspace(-10, 10, 16))

		xx.setflags(write = False)
		yy.setflags(write = False)
		zz.setflags(write = False)

		potential = electrostatic(graphite, xx, yy, zz)

class TestPElectrostatic(unittest.TestCase):

	def setUp(self):
		self.crystal = deepcopy(graphite)

	def test_return_shape(self):
		""" Test that the return shape of pelectrostatic is the same as input arrays """
		xx, yy = np.meshgrid(np.linspace(-10, 10, 32), np.linspace(-10, 10, 32))
		potential = pelectrostatic(graphite, xx, yy)

		self.assertSequenceEqual(xx.shape, potential.shape)
	
	def test_side_effects(self):
		""" Test that mesh arrays are not written to in pelectrostatic """
		xx, yy = np.meshgrid(np.linspace(-10, 10, 32), np.linspace(-10, 10, 32))
		xx.setflags(write = False)
		yy.setflags(write = False)
		potential = pelectrostatic(graphite, xx, yy)
	
	def test_trivial(self):
		""" Test that the projected electrostatic potential from an empty slice of crystal is zero"""
		xx, yy = np.meshgrid(np.linspace(-10, 10, 32), np.linspace(-10, 10, 32))
		potential = pelectrostatic(graphite, xx, yy, bounds = (1,2))
		self.assertTrue(np.allclose(potential, 0))

if __name__ == '__main__':
	unittest.main()