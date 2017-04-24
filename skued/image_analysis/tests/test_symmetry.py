
import numpy as np
from .. import angular_average
import unittest

class TestAngularAverage(unittest.TestCase):

	def test_trivial_array(self):
		""" Test angular_average on an array of zeroes """
		image = np.zeros(shape = (256, 256), dtype = np.float)
		center = (image.shape[0]/2, image.shape[1]/2)

		extras = dict()
		intensity = angular_average(image, center, extras = extras)

		self.assertTrue(intensity.sum() == 0)

		radius, error = extras['radius'], extras['error']
		self.assertSequenceEqual(radius.shape, error.shape)
		self.assertSequenceEqual(intensity.shape, radius.shape)
		
	def test_ring(self):
		""" Test angular_average on an image with a wide ring """
		image = np.zeros(shape = (256, 256), dtype = np.float)
		center = (image.shape[0]/2, image.shape[1]/2)
		xc, yc = center

		# Create an image with a wide ring
		extent = np.arange(0, image.shape[0])
		xx, yy = np.meshgrid(extent, extent)
		rr = np.sqrt((xx - xc)**2 + (yy - yc)**2)
		image[np.logical_and(24 < rr,rr < 26)] = 1

		intensity = angular_average(image, center)
		self.assertEqual(intensity.max(), image.max())

if __name__ == '__main__':
	unittest.main()