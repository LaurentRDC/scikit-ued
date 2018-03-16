# -*- coding: utf-8 -*-

import numpy as np
from .. import azimuthal_average, powder_center
from ... import Crystal, powdersim
import unittest

from skimage.filters import gaussian

np.random.seed(23)

def circle_image(shape, center, radii, intensities):
    """ Creates an image with circle or thickness 2 """
    im = np.zeros(shape = shape, dtype = np.float)
    xx, yy = np.ogrid[0:shape[0],0:shape[1]]
    xx, yy = xx - center[0], yy - center[1]
    for radius, intensity in zip(radii, intensities):
        rr = np.sqrt(xx**2 + yy**2)
        im[np.logical_and(rr < radius + 1, rr > radius - 1)] = intensity
    
    im[:] = gaussian(im, 5)
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
        self.assertTrue(np.allclose(center, powder_center(im), atol = 1))

class TestAzimuthalAverage(unittest.TestCase):
    
    def test_trivial_array(self):
        """ Test azimuthal_average on an array of zeroes """
        image = np.zeros(shape = (256, 256), dtype = np.float)
        center = (image.shape[0]/2, image.shape[1]/2)

        radius, intensity = azimuthal_average(image, center)
        
        self.assertTrue(intensity.sum() == 0)
        self.assertSequenceEqual(intensity.shape, radius.shape)
        
    def test_ring(self):
        """ Test azimuthal_average on an image with a wide ring """
        image = np.zeros(shape = (256, 256), dtype = np.float)
        center = (image.shape[0]/2, image.shape[1]/2)
        xc, yc = center

        # Create an image with a wide ring
        extent = np.arange(0, image.shape[0])
        xx, yy = np.meshgrid(extent, extent)
        rr = np.sqrt((xx - xc)**2 + (yy - yc)**2)
        image[np.logical_and(24 < rr,rr < 26)] = 1

        radius, intensity = azimuthal_average(image, center)
        self.assertEqual(intensity.max(), image.max())
        self.assertSequenceEqual(radius.shape, intensity.shape)
    
    def test_angular_bounds(self):
        """ Test azimuthal_average with a restrictive angular_bounds argument """
        image = np.zeros(shape = (256, 256), dtype = np.float)
        center = (image.shape[0]/2, image.shape[1]/2)
        xc, yc = center

        # Create an image with a wide ring
        extent = np.arange(0, image.shape[0])
        xx, yy = np.meshgrid(extent, extent)
        rr = np.sqrt((xx - xc)**2 + (yy - yc)**2)
        angles = np.rad2deg(np.arctan2(yy - yc, xx - xc)) + 180
        image[np.logical_and(0 <= angles, angles <= 60)] = 1

        with self.subTest('0 - 360'):
            radius, intensity = azimuthal_average(image, center, angular_bounds = None)
            r360, int360 = azimuthal_average(image, center, angular_bounds = (0, 360))
            self.assertTrue(np.allclose(intensity, int360))
        
        with self.subTest('Inside angle bounds'):
            radius, intensity = azimuthal_average(image, center, angular_bounds = (0, 60))
            self.assertTrue(np.allclose(intensity, np.ones_like(intensity)))

        with self.subTest('Overlapping bounds'):
            radius, intensity = azimuthal_average(image, center, angular_bounds = (15, 75))
            self.assertFalse(np.all(intensity < np.ones_like(intensity)))
        
        with self.subTest('Outside angle bounds'):
            radius, intensity = azimuthal_average(image, center, angular_bounds = (60, 360))
            self.assertTrue(np.allclose(intensity, np.zeros_like(intensity)))
        
        with self.subTest('Inside angle bounds with 360deg rollover'):
            radius, intensity = azimuthal_average(image, center, angular_bounds = (60 + 360, 360 + 360))
            self.assertTrue(np.allclose(intensity, np.zeros_like(intensity)))
    
    def test_ring_with_mask(self):
        """ Test azimuthal_average on an image with a wide ring """
        image = np.zeros(shape = (256, 256), dtype = np.float)
        center = (image.shape[0]/2, image.shape[1]/2)
        xc, yc = center

        mask = np.zeros_like(image, dtype = np.bool)
        mask[120:140, 0:140] = True

        # Create an image with a wide ring
        extent = np.arange(0, image.shape[0])
        xx, yy = np.meshgrid(extent, extent)
        rr = np.sqrt((xx - xc)**2 + (yy - yc)**2)
        image[np.logical_and(24 < rr,rr < 26)] = 1

        radius, intensity = azimuthal_average(image, center, mask = mask)

        self.assertEqual(intensity.max(), image.max())
        self.assertSequenceEqual(radius.shape, intensity.shape)
    
    def test_trim_and_mask(self):
        """ Test that regions that only have masks contributions are not present
        in the angular average """
        image = np.ones(shape = (256, 256), dtype = np.float)
        center = (image.shape[0]/2, image.shape[1]/2)
        xc, yc = center

        # Create an image with a wide ring
        extent = np.arange(0, image.shape[0])
        xx, yy = np.meshgrid(extent, extent)
        rr = np.hypot(xx - xc, yy - yc)

        mask = np.zeros_like(image, dtype = np.bool)
        mask[rr < 20] = True
        #image[rr < 20] = 0

        radius, intensity = azimuthal_average(image, center, mask = mask, trim = False)
        self.assertEqual(radius.min(), 0)

        radius_trimmed, intensity_trimmed = azimuthal_average(image, center, mask = mask, trim = True)
        self.assertEqual(radius_trimmed.min(), 20)
    
    def test_mask_and_nan(self):
        """ Test that azimuthal_average with masks does not yield NaNs. This can happen for large masks. """
        image = np.ones(shape = (256, 256), dtype = np.int16)
        mask = np.zeros_like(image, dtype = np.bool)
        mask[100:156, 100:156] = True
        _, av = azimuthal_average(image, center = (128, 128), mask = mask, trim = False)

        self.assertFalse(np.any(np.isnan(av)))


if __name__ == '__main__':
    unittest.main()