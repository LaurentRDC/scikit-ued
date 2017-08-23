# -*- coding: utf-8 -*-
import numpy as np
from .. import mask_from_collection, combine_masks, mask_image
import unittest

class TestMaskFromCollection(unittest.TestCase):

    def test_trivial(self):
        """ Test on set of images with value zero """
        images = [np.zeros((64,64)) for _ in range(5)]
        mask = mask_from_collection(images)

        self.assertSequenceEqual(images[0].shape, mask.shape)
        self.assertFalse(np.any(mask))

    def test_intensity_threshold_no_lower_bound(self):
        """ Test that intensity threshold is respected """
        images = [np.zeros((64,64)) for _ in range(5)]
        images[2][32,4] = 10
        mask = mask_from_collection(images, px_thresh = 9)

        self.assertEqual(np.sum(mask), 1)   # only one pixels is masked
        self.assertEqual(mask[32, 4], True)

    def test_intensity_threshold_with_lower_bound(self):
        """ Test that intensity threshold is respected """
        images = [np.zeros((64,64)) for _ in range(5)]
        images[2][32,4] = -10
        mask = mask_from_collection(images, px_thresh = (0, np.inf))

        self.assertEqual(np.sum(mask), 1)   # only one pixels is masked
        self.assertEqual(mask[32, 4], True)
    
    def test_std_threshold(self):
        """ Test that std threshold is respected """
        images = [np.zeros((64,64)) for _ in range(5)]
        images[0][5, 12] = 1000
        mask = mask_from_collection(images, px_thresh = 10000, std_thresh = 1)

        self.assertEqual(np.sum(mask), 1)   # only one pixels is masked
        self.assertEqual(mask[5, 12], True)
    
    def test_single_image(self):
        """ Test that mask_from_collection works even if input
        is a single array """
        images = np.zeros((64,64))
        mask = mask_from_collection(images)

        self.assertFalse(np.any(mask))

class TestCombineMasks(unittest.TestCase):

    def test_trivial(self):
        masks = tuple([np.zeros((64,64)) for _ in range(5)])
        combined = combine_masks(*masks)

        self.assertFalse(np.any(combined))
    
    def test_single_element(self):
        masks = tuple([np.zeros((64,64)) for _ in range(5)])
        masks[0][4, 6] = True
        combined = combine_masks(*masks)

        self.assertEqual(np.sum(combined), 1)
        self.assertTrue(combined[4, 6])

class TestMaskImage(unittest.TestCase):

    def test_trivial(self):
        mask = np.zeros((64,64), dtype = np.bool)
        image = np.random.random((64,64))
        masked = mask_image(image, mask)

        self.assertTrue(np.allclose(image, masked))
    
    def test_no_copy(self):
        """ Test that mask_image can work in-place """
        mask = np.random.randint(0, 1, size = (64,64), dtype = np.bool)
        image = np.random.random((64,64))
        masked = mask_image(image, mask, copy = False)

        self.assertIs(image, masked)

if __name__ == '__main__':
    unittest.main()