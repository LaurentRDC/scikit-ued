# -*- coding: utf-8 -*-
import unittest
from pathlib import Path
from random import randint

import numpy as np
from skimage import data
from skimage.filters import gaussian
from skimage.io import imread
from skimage.feature import register_translation
from skimage.transform import rotate
from scipy.ndimage import fourier_shift

from .. import (
    align,
    ialign,
    itrack_peak,
    masked_register_translation,
    shift_image,
)
from .test_powder import circle_image

np.random.seed(23)


class TestIAlign(unittest.TestCase):
    def test_trivial(self):
        """ Test alignment of identical images """
        aligned = tuple(ialign([data.camera() for _ in range(5)]))

        self.assertEqual(len(aligned), 5)
        self.assertSequenceEqual(data.camera().shape, aligned[0].shape)

    def test_misaligned_canned_images_fast(self):
        """ shift images from skimage.data by entire pixels.
	   We don't expect perfect alignment."""
        original = data.camera()
        misaligned = [
            shift_image(original, (randint(-4, 4), randint(-4, 4))) for _ in range(5)
        ]

        aligned = ialign(misaligned, reference=original, fast=True)

        # TODO: find a better figure-of-merit for alignment
        for im in aligned:
            # edge will be filled with zeros, we ignore
            diff = np.abs(original[5:-5, 5:-5] - im[5:-5, 5:-5])

            # Want less than 1% difference
            percent_diff = np.sum(diff) / (
                diff.size * (original.max() - original.min())
            )
            self.assertLess(percent_diff, 1)

    def test_misaligned_canned_images_notfast(self):
        """ shift images from skimage.data by entire pixels.
	   We don't expect perfect alignment."""
        original = data.camera()
        misaligned = [
            shift_image(original, (randint(-4, 4), randint(-4, 4))) for _ in range(5)
        ]

        aligned = ialign(misaligned, reference=original, fast=False)

        # TODO: find a better figure-of-merit for alignment
        for im in aligned:
            # edge will be filled with zeros, we ignore
            diff = np.abs(original[5:-5, 5:-5] - im[5:-5, 5:-5])

            # Want less than 1% difference
            percent_diff = np.sum(diff) / (
                diff.size * (original.max() - original.min())
            )
            self.assertLess(percent_diff, 1)


class TestShiftImage(unittest.TestCase):
    def test_trivial(self):
        """ Shift an array by (0,0) """
        arr = np.random.random(size=(64, 64))
        shifted = shift_image(arr, (0, 0))
        self.assertTrue(np.allclose(arr, shifted))

    def test_shift_float16(self):
        """ Interpolation requires float32 or more bits. """
        arr = np.random.random(size=(64, 64)).astype(np.float16)
        shifted1 = shift_image(arr, (5, -3))
        shifted2 = shift_image(shifted1, (-5, 3))
        self.assertTrue(np.allclose(arr[5:-5, 5:-5], shifted2[5:-5, 5:-5]))

    def test_back_and_forth(self):
        """ Test shift_image in two directions """
        arr = np.random.random(size=(64, 64))
        shifted1 = shift_image(arr, (5, -3))
        shifted2 = shift_image(shifted1, (-5, 3))
        self.assertTrue(np.allclose(arr[5:-5, 5:-5], shifted2[5:-5, 5:-5]))

    def test_return_type(self):
        """ Test that a shifted array will cast accordingly to the fill_value """
        arr = np.random.randint(0, 255, size=(64, 64), dtype=np.uint8)
        shifted = shift_image(arr, shift=(10, 10), fill_value=np.nan)  # np.nan is float
        self.assertEqual(shifted.dtype, np.float)

    def test_out_of_bounds(self):
        """ Test that shifting by more than the size of an array
		returns an array full of the fill_value parameter """
        arr = np.random.random(size=(64, 64)) + 1  # no zeros in this array
        shifted = shift_image(arr, (128, 128), fill_value=0.0)
        self.assertTrue(np.allclose(np.zeros_like(arr), shifted))

    def test_fill_value(self):
        """ Test that shifted array edges are filled with the correct value """
        arr = np.random.random(size=(64, 64))
        shifted = shift_image(arr, shift=(0, 10), fill_value=np.nan)
        self.assertTrue(np.all(np.isnan(shifted[:10, :])))


class TestAlign(unittest.TestCase):
    def test_no_side_effects(self):
        """ Test that aligned images are not modified in-place """
        im = np.array(data.camera()[0:64, 0:64])
        im.setflags(write=False)
        aligned = align(im, reference=im, fill_value=np.nan)
        self.assertEqual(im.dtype, data.camera().dtype)

    def test_misaligned_canned_images_fast(self):
        """ shift images from skimage.data by entire pixels.
	   	We don't expect perfect alignment."""
        original = data.camera()
        misaligned = shift_image(original, (randint(-4, 4), randint(-4, 4)))

        aligned = align(misaligned, reference=original, fast=True)

        # edge will be filled with zeros, we ignore
        diff = np.abs(original[5:-5, 5:-5] - aligned[5:-5, 5:-5])

        # Want less than 1% difference
        percent_diff = np.sum(diff) / (diff.size * (original.max() - original.min()))
        self.assertLess(percent_diff, 1)

    def test_misaligned_canned_images_notfast(self):
        """ shift images from skimage.data by entire pixels.
	   	We don't expect perfect alignment."""
        original = data.camera()
        misaligned = shift_image(original, (randint(-4, 4), randint(-4, 4)))

        aligned = align(misaligned, reference=original, fast=False)

        # edge will be filled with zeros, we ignore
        diff = np.abs(original[5:-5, 5:-5] - aligned[5:-5, 5:-5])

        # Want less than 1% difference
        percent_diff = np.sum(diff) / (diff.size * (original.max() - original.min()))
        self.assertLess(percent_diff, 1)


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


class TestMaskedRegisterTranslation(unittest.TestCase):
    def test_padfield_data_full_mode(self):
        """ Test translation registration for data included in Padfield 2010 in `full` convolution mode"""
        # Test translated from MATLABimplementation `MaskedFFTRegistrationTest` file. You can find the source code here:
        # http://www.dirkpadfield.com/Home/MaskedFFTRegistrationCode.zip
        IMAGES_DIR = Path(__file__).parent / "images"

        shifts = [(75, 75), (-130, 130), (130, 130)]
        for xi, yi in shifts:
            with self.subTest("X = {:d}, Y = {:d}".format(xi, yi)):
                fixed_image = imread(
                    IMAGES_DIR / "OriginalX{:d}Y{:d}.png".format(xi, yi)
                )
                moving_image = imread(
                    IMAGES_DIR / "TransformedX{:d}Y{:d}.png".format(xi, yi)
                )

                # Valid pixels are 1
                fixed_mask = fixed_image != 0
                moving_mask = moving_image != 0

                # Note that shifts in x and y and shifts in cols and rows
                shift_y, shift_x = masked_register_translation(
                    fixed_image,
                    moving_image,
                    fixed_mask,
                    moving_mask,
                    mode="full",
                    overlap_ratio=1 / 10,
                )
                # NOTE: by looking at the test code from Padfield's MaskedFFTRegistrationCode repository,
                # 		the shifts were not xi and yi, but xi and -yi
                self.assertTupleEqual((-xi, yi), (shift_x, shift_y))

    def test_padfield_data_same_mode(self):
        """ Test translation registration for data included in Padfield 2010 in `same` convolution mode"""
        # Test translated from MATLABimplementation `MaskedFFTRegistrationTest` file. You can find the source code here:
        # http://www.dirkpadfield.com/Home/MaskedFFTRegistrationCode.zip
        IMAGES_DIR = Path(__file__).parent / "images"

        shifts = [(75, 75), (-130, 130), (130, 130)]
        for xi, yi in shifts:
            with self.subTest("X = {:d}, Y = {:d}".format(xi, yi)):
                fixed_image = imread(
                    IMAGES_DIR / "OriginalX{:d}Y{:d}.png".format(xi, yi)
                )
                moving_image = imread(
                    IMAGES_DIR / "TransformedX{:d}Y{:d}.png".format(xi, yi)
                )

                # Valid pixels are 1
                fixed_mask = fixed_image != 0
                moving_mask = moving_image != 0

                # Note that shifts in x and y and shifts in cols and rows
                shift_y, shift_x = masked_register_translation(
                    fixed_image,
                    moving_image,
                    fixed_mask,
                    moving_mask,
                    mode="same",
                    overlap_ratio=1 / 10,
                )
                # NOTE: by looking at the test code from Padfield's MaskedFFTRegistrationCode repository,
                # 		the shifts were not xi and yi, but xi and -yi
                self.assertTupleEqual((-xi, yi), (shift_x, shift_y))

    def test_against_skimage_register_translation(self):
        """ Test masked_register_translation against scikit-image's register_translation for trivial masks """

        # Generate shifted image like scikit-image' test suite
        # Most importantly, image shifting is done the same way, using a Fourier shift filter.
        # https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/tests/test_register_translation.py
        reference_image = data.camera()
        shift = (-7, 12)
        shifted = np.real(
            np.fft.ifft2(fourier_shift(np.fft.fft2(reference_image), shift))
        )
        trivial_mask = np.ones_like(reference_image, dtype=np.bool)

        nonmasked_result, *_ = register_translation(reference_image, shifted)
        masked_result = masked_register_translation(
            reference_image, shifted, trivial_mask, overlap_ratio=1 / 10
        )

        self.assertTrue(np.allclose(nonmasked_result, masked_result))


if __name__ == "__main__":
    unittest.main()
