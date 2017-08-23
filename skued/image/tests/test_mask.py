# -*- coding: utf-8 -*-
import numpy as np
from .. import mask_from_collection
import unittest

class TestMaskFromCollection(unittest.TestCase):

    def test_trivial(self):
        """ Test on set of images with value zero """
        images = [np.zeros((64,64)) for _ in range(5)]
        mask = mask_from_collection(images)

        self.assertSequenceEqual(images[0].shape, mask.shape)
        self.assertFalse(np.any(mask))

    def test_intensity_threshold(self):
        """ Test that intensity threshold is respected """
        images = [np.zeros((64,64)) for _ in range(5)]
        images[2][32,4] = 10
        mask = mask_from_collection(images, int_thresh = 9)

        self.assertEqual(np.sum(mask), 1)   # only one pixels is masked
        self.assertEqual(mask[32, 4], True)
    
    def test_std_threshold(self):
        """ Test that std threshold is respected """
        images = [np.zeros((64,64)) for _ in range(5)]
        images[0][5, 12] = 1000
        mask = mask_from_collection(images, int_thresh = 10000, std_thresh = 1)

        self.assertEqual(np.sum(mask), 1)   # only one pixels is masked
        self.assertEqual(mask[5, 12], True)
    
    def test_single_image(self):
        """ Test that mask_from_collection works even if input
        is a single array """
        images = np.zeros((64,64))
        mask = mask_from_collection(images)

        self.assertFalse(np.any(mask))

if __name__ == '__main__':
    unittest.main()