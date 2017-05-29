# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from .. import angular_average, powder_center
import unittest

from skimage.filters import gaussian

np.random.seed(23)

def circle_image(shape, center, radii, intensities):
	""" Creates an image with circle or thickness 2 """
	im = np.zeros(shape = shape, dtype = np.float)
	xx, yy = np.meshgrid(np.arange(0, shape[0]),
						 np.arange(0, shape[1]),
						 indexing = 'ij')
	xx, yy = xx - center[0], yy - center[1]
	for radius, intensity in zip(radii, intensities):
		rr = np.sqrt(xx**2 + yy**2)
		im[np.logical_and(rr < radius + 1, rr > radius - 1)] = intensity
	
	return im

class TestPowderCenter(unittest.TestCase):
	
	def test_trivial(self):
		""" Test center-finding without any noise """
		center = (64, 64)
		im = circle_image(shape = (128, 128), center = center, 
						  radii = [16, 32], intensities = [2,1])
		self.assertSequenceEqual(center, powder_center(im))

	def test_with_noise(self):
		""" Test center-finding with noise """
		center = (58, 67)
		im = circle_image(shape = (128, 128), center = center, 
						  radii = [16, 32], intensities = [4,3])
		im += (im.max() / 5) * np.random.random(size = im.shape)
		self.assertSequenceEqual(center, powder_center(im))
	
	def test_with_low_signal(self):
		center = (250, 300)
		im = circle_image(shape = (512, 512), center = center, 
						  radii = [32, 125, 200],
						  intensities = [40, 30, 40])
		im += (im.max() / 10) * np.random.random(size = im.shape)
		im[:] = gaussian(im, sigma = 5)
		self.assertSequenceEqual(center, powder_center(im))

class TestAngularAverage(unittest.TestCase):

	def test_trivial_array(self):
		""" Test angular_average on an array of zeroes """
		image = np.zeros(shape = (256, 256), dtype = np.float)
		center = (image.shape[0]/2, image.shape[1]/2)

		extras = dict()
		radius, intensity = angular_average(image, center, extras = extras)
		
		self.assertTrue(intensity.sum() == 0)

		error = extras['error']
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

		radius, intensity = angular_average(image, center)
		self.assertEqual(intensity.max(), image.max())
	
	def test_ring_with_mask(self):
		""" Test angular_average on an image with a wide ring """
		image = np.zeros(shape = (256, 256), dtype = np.float)
		center = (image.shape[0]/2, image.shape[1]/2)
		xc, yc = center

		mask = np.zeros_like(image, dtype = np.bool)
		mask[120:140, 0:128] = True

		# Create an image with a wide ring
		extent = np.arange(0, image.shape[0])
		xx, yy = np.meshgrid(extent, extent)
		rr = np.sqrt((xx - xc)**2 + (yy - yc)**2)
		image[np.logical_and(24 < rr,rr < 26)] = 1

		radius, intensity = angular_average(image, center, mask = mask)

		self.assertEqual(intensity.max(), image.max())

if __name__ == '__main__':
	unittest.main()