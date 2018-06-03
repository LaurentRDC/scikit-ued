# -*- coding: utf-8 -*-
import numpy as np
from .. import nfold, reflection
import unittest
from warnings import catch_warnings, simplefilter

np.random.seed(23)

class TestNFoldSymmetry(unittest.TestCase):

    def test_trivial(self):
        """ Test nfold_symmetry averaging on trivial array """
        im = np.zeros( (256, 256) )
        rot = nfold(im, mod = 3)
        self.assertTrue(np.allclose(rot, im))

    def test_valid_mod(self):
        """ Test the the N-fold symmetry argument is valid """
        im = np.empty( (128, 128) )
        with self.assertRaises(ValueError):
            nfold(im, mod = 1.7)
    
    def test_mask(self):
        """ Test that nfold_symmetry() works correctly with a mask """
        im = np.zeros((128, 128), dtype = np.int)
        mask = np.zeros_like(im, dtype = np.bool)

        with catch_warnings():
            simplefilter('ignore')
            im[0:20] = 1
            mask[0:20] = True
            
            rot = nfold(im, mod = 2, mask = mask)
            self.assertTrue(np.allclose(rot, np.zeros_like(rot)))
    
    def test_no_side_effects(self):
        """ Test that nfold() does not modify the input image and mask """
        im = np.empty((128, 128), dtype = np.float)
        mask = np.zeros_like(im, dtype = np.bool)

        im.setflags(write = False)
        mask.setflags(write = False)

        rot = nfold(im, center = (67, 93),mod = 3, mask = mask)
    
    def test_fill_value(self):
        """ Test that the fill_value parameter of nfold() is working correctly """
        im = 1000*np.random.random(size = (256, 256))
        mask = np.random.choice([True, False], size = im.shape)

        with self.subTest('fill_value = np.nan'):
            with catch_warnings():
                simplefilter('ignore')
                rot = nfold(im, center = (100, 150), mod = 5, mask = mask, fill_value = np.nan)

            self.assertTrue(np.any(np.isnan(rot)))

        with self.subTest('fill_value = 0.0'):
            with catch_warnings():
                simplefilter('ignore')
                rot = nfold(im, center = (100, 150), mod = 5, mask = mask, fill_value = 0.0)

            self.assertFalse(np.any(np.isnan(rot)))
    
    def test_output_range(self):
        """ Test that nfold() does not modify the value range """
        im = 1000*np.random.random(size = (256, 256))
        mask = np.random.choice([True, False], size = im.shape)

        with catch_warnings():
            simplefilter('ignore')
            rot = nfold(im, center = (100, 150), mod = 5, mask = mask)

        self.assertLessEqual(rot.max(), im.max())
        # In the case of a mask that overlaps with itself when rotated,
        # the average will be zero due to nan_to_num
        self.assertGreaterEqual(rot.min(), min(im.min(), 0))
    
    def test_mod_1(self):
        """ Test that nfold(mod = 1) returns an unchanged image, except
        perhaps for a cast to float """
        im = 1000*np.random.random(size = (256, 256))
        rot = nfold(im, mod = 1)
        self.assertTrue(np.allclose(im, rot))

class TestReflectionSymmetry(unittest.TestCase):

    def test_trivial(self):
        """ Test the reflection symmetry of an image of zeroes """
        im = np.zeros( (256, 256) )
        ref = reflection(im, angle = -15.3)
        self.assertTrue(np.allclose(ref, im))

    def test_no_side_effects(self):
        """ Test that reflection() does not modify the input image and mask """
        im = np.empty((128, 128), dtype = np.float)
        mask = np.zeros_like(im, dtype = np.bool)

        im.setflags(write = False)
        mask.setflags(write = False)

        rot = reflection(im, center = (67, 93),angle = 35, mask = mask)

    def test_output_range(self):
        """ Test that reflection() does not modify the value range """
        im = 1000*np.random.random(size = (256, 256))
        mask = np.random.choice([True, False], size = im.shape)

        with catch_warnings():
            simplefilter('ignore')
            rot = reflection(im, center = (100, 150),angle = 5, mask = mask)

        self.assertLessEqual(rot.max(), im.max())
        # In the case of a mask that overlaps with itself when rotated,
        # the average will be zero due to nan_to_num
        self.assertGreaterEqual(rot.min(), min(im.min(), 0))
    
    def test_correctness_angle0(self):
        """ Test that reflection() correctly symmetrizes around the x-axis. """

        im = np.zeros((256, 256), dtype = np.float)
        im[0:10,:] = 1

        reflected = reflection(im, angle = 0)

        expected = np.array(im, copy = True)
        expected[0:10, :] = 0.5
        expected[246:256, :] = 0.5

        self.assertTrue(np.allclose(expected, reflected))
    
    def test_angle_vs_angle_plus_180(self):
        """ Test that the result of reflection() is the same for any angle and angle + 180 """
        im = np.random.random((256,256))
        reflected1 = reflection(im, angle = 15)
        reflected2 = reflection(im, angle = 195)    # + 180

        self.assertTrue(np.allclose(reflected1, reflected2))

    def test_correctness_angle90(self):
        """ Test that reflection() correctly symmetrizes around the y-axis. """

        im = np.zeros((256, 256), dtype = np.float)
        im[:, 0:10] = 1

        reflected = reflection(im, angle = 90)

        expected = np.array(im, copy = True)
        expected[:, 0:10] = 0.5
        expected[:, 246:256] = 0.5

        self.assertTrue(np.allclose(expected, reflected))

if __name__ == '__main__':
    unittest.main()