# -*- coding: utf-8 -*-
import unittest
from pathlib import Path
from random import randint

import numpy as np
from scipy import ndimage as ndi
from skimage import data
from skimage.feature import register_translation
from skimage.filters import gaussian
from skimage.transform import rotate
from skimage.data import camera

from skued import align, ialign, itrack_peak

from .test_powder import circle_image

np.random.seed(23)

TEST_IMAGE = camera()[::2, ::2]  # decimate by 2 for faster tests


class TestIAlign(unittest.TestCase):
    def test_trivial(self):
        """ Test alignment of identical images """
        aligned = tuple(ialign([TEST_IMAGE for _ in range(5)]))

        self.assertEqual(len(aligned), 5)
        self.assertSequenceEqual(TEST_IMAGE.shape, aligned[0].shape)

    def test_misaligned_canned_images(self):
        """shift images from skimage.data by entire pixels.
        We don't expect perfect alignment."""
        original = TEST_IMAGE
        misaligned = [original] + [
            ndi.shift(original, (randint(-4, 4), randint(-4, 4))) for _ in range(5)
        ]
        # Mask the edges
        mask = np.ones_like(TEST_IMAGE, dtype=np.bool)
        mask[0:5, :] = False
        mask[:-5, :] = False
        mask[:, 0:5] = False
        mask[:, :-5] = False

        aligned = ialign(misaligned, reference=original, mask=mask)

        # TODO: find a better figure-of-merit for alignment
        for im in aligned:
            # edge will be filled with zeros, we ignore
            diff = np.abs(original[mask] - im[mask])

            # Want less than 0.5% difference
            percent_diff = np.sum(diff) / (
                diff.size * (original.max() - original.min())
            )
            self.assertLess(percent_diff, 1)


class TestAlign(unittest.TestCase):
    def test_no_side_effects(self):
        """ Test that aligned images are not modified in-place """
        im = np.array(TEST_IMAGE[0:64, 0:64])
        im.setflags(write=False)
        aligned = align(im, reference=im, fill_value=np.nan)
        self.assertEqual(im.dtype, TEST_IMAGE.dtype)

    def test_misaligned_canned_images(self):
        """shift images from skimage.data by entire pixels.
        We don't expect perfect alignment."""
        original = TEST_IMAGE
        misaligned = ndi.shift(original, (-4, 4))

        # Mask the edges
        mask = np.ones_like(TEST_IMAGE, dtype=np.bool)
        mask[0:5, :] = False
        mask[:-5, :] = False
        mask[:, 0:5] = False
        mask[:, :-5] = False

        aligned = align(misaligned, reference=original, mask=mask)

        # edge will be filled with zeros, we ignore
        diff = np.abs(original[mask] - aligned[mask])

        # Want less than 0.5% difference
        percent_diff = np.sum(diff) / (diff.size * (original.max() - original.min()))
        self.assertLess(percent_diff, 0.5)


class TestItrackPeak(unittest.TestCase):
    def test_trivial(self):
        """ Test that shift is identically zero for images that are identical """
        # Array prototype is just zeros
        # with a 'peak' in the center
        prototype = np.zeros(shape=(17, 17))
        prototype[9, 9] = 10
        images = [np.array(prototype) for _ in range(20)]
        shifts = itrack_peak(images, row_slice=np.s_[:], col_slice=np.s_[:])

        for shift in shifts:
            self.assertTrue(np.allclose(shift, (0.0, 0.0)))

    def test_length(self):
        """ Test that shifts yielded by itrack_peak are as numerous as the number of input pictures """
        images = [np.random.random(size=(4, 4)) for _ in range(20)]
        shifts = list(itrack_peak(images, row_slice=np.s_[:], col_slice=np.s_[:]))

        self.assertEqual(len(shifts), len(images))


if __name__ == "__main__":
    unittest.main()
