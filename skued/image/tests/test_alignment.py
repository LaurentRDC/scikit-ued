# -*- coding: utf-8 -*-
import unittest
from pathlib import Path
from random import randint

import numpy as np
from scipy import ndimage as ndi
from skimage import data
from skimage.draw import ellipse
from skimage.feature import register_translation
from skimage.filters import gaussian
from skimage.transform import rotate
import skimage.data as data

from skued import align, ialign, itrack_peak

from .test_powder import circle_image

np.random.seed(23)


def camera():
    return data.camera()[::2, ::2]  # decimate by 2 for faster tests


class TestIAlign(unittest.TestCase):
    def test_trivial(self):
        """ Test alignment of identical images """
        aligned = tuple(ialign([camera() for _ in range(5)]))

        self.assertEqual(len(aligned), 5)
        self.assertSequenceEqual(camera().shape, aligned[0].shape)

    def test_misaligned_canned_images(self):
        """shift images from skimage.data by entire pixels.
        We don't expect perfect alignment."""
        reference = camera()
        shifts = [(randint(-4, 0), randint(0, 4)) for _ in range(5)]

        mask = np.zeros_like(reference, dtype=np.bool)
        rr1, cc1 = ellipse(129, 127, r_radius=63, c_radius=50, shape=reference.shape)
        mask[rr1, cc1] = True

        misaligned = (ndi.shift(camera(), shift=s) for s in shifts)

        for aligned, (sx, sy) in zip(
            ialign(misaligned, reference=reference, mask=mask), shifts
        ):
            self.assertTrue(np.allclose(reference[sx::, 0:-sy], aligned[sx::, 0:-sy]))


class TestAlign(unittest.TestCase):
    def test_no_side_effects(self):
        """ Test that aligned images are not modified in-place """
        im = np.array(camera()[0:64, 0:64])
        im.setflags(write=False)
        aligned = align(im, reference=im, fill_value=np.nan)
        self.assertEqual(im.dtype, camera().dtype)

    def test_misaligned_no_mask(self):
        """ Test that alignment of images with no masks works """
        reference = camera()
        image = ndi.shift(camera(), shift=(-7, 12))
        aligned = align(image, reference=reference)

        self.assertTrue(np.allclose(reference[7::, 0:-12], aligned[7::, 0:-12]))

    def test_misaligned_with_mask(self):
        """ Test that alignment of images with no masks works """
        reference = camera()
        image = ndi.shift(camera(), shift=(-7, 12))

        mask = np.ones_like(reference, dtype=np.bool)
        mask[75:100, 50:100] = False

        aligned = align(image, reference=reference, mask=mask)

        self.assertTrue(np.allclose(reference[7::, 0:-12], aligned[7::, 0:-12]))


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
