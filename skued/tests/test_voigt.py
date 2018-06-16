# -*- coding: utf-8 -*-
import numpy as np
import unittest
from .. import gaussian, lorentzian, pseudo_voigt

class TestFunctionUnitIntegral(object):
	# Same tests for lorentzian and gaussian
	# Unit integral testing is only to one decimal place
	# to save some time.

	def test_1d_unit_integral(self):
		dx = 0.1
		x = np.arange(-5,5,dx)
		center = 0

		integral = np.sum(self.function(x, center, 0.2))*dx
		self.assertAlmostEqual(integral, 1, places = 1)

	def test_2d_unit_integral(self):
		dx = dy = 0.1
		extent = np.arange(-5,5,dx)
		xx, yy = np.meshgrid(extent,extent)
		center = [0,0]

		integral = np.sum(self.function([xx, yy], center, 0.2))*dx*dy
		self.assertAlmostEqual(integral, 1, places = 1)

	def test_3d_unit_integral(self):
		dx = dy = dz = 0.1
		extent = np.arange(-3,3,dx)
		xx, yy, zz = np.meshgrid(extent,extent,extent)
		center = [0,0, 0]

		integral = np.sum(self.function([xx, yy, zz], center, 0.2))*dx*dy*dz
		self.assertAlmostEqual(integral, 1, places = 1)

class TestGaussianUnitIntegral(TestFunctionUnitIntegral, unittest.TestCase):
	def setUp(self):
		 self.function = gaussian

class TestLorentzianUnitIntegral(TestFunctionUnitIntegral, unittest.TestCase):
	def setUp(self):
		 self.function = lorentzian

class TestVoigtFunctionUnitIntegral(unittest.TestCase):
	# voigt and pseudo_voigt have two width parameters
	# Unit integral testing is only to one decimal place
	# to save some time.

	def setUp(self):
		self.function = pseudo_voigt

	def test_1d_unit_integral(self):
		dx = 0.1
		x = np.arange(-5,5,dx)
		center = 0

		integral = np.sum(self.function(x, center, 0.2, 0.3))*dx
		self.assertAlmostEqual(integral, 1, places = 1)

	def test_2d_unit_integral(self):
		dx = dy = 0.1
		extent = np.arange(-5,5,dx)
		xx, yy = np.meshgrid(extent,extent)
		center = [0,0]

		integral = np.sum(self.function([xx, yy], center, 0.2, 0.3))*dx*dy
		self.assertAlmostEqual(integral, 1, places = 1)

	def test_3d_unit_integral(self):
		dx = dy = dz = 0.1
		extent = np.arange(-3,3,dx)
		xx, yy, zz = np.meshgrid(extent,extent,extent)
		center = [0,0, 0]

		integral = np.sum(self.function([xx, yy, zz], center, 0.2, 0.3))*dx*dy*dz
		self.assertAlmostEqual(integral, 1, places = 1)

if __name__ == '__main__':
	unittest.main()