# -*- coding: utf-8 -*-
import unittest
from pathlib import Path

import numpy as np
from scipy.signal import correlate
from skimage.io import imread
from skimage.data import camera

from .. import mnxc, xcorr

np.random.seed(23)


class TestXcorr(unittest.TestCase):
    def test_final_shape(self):
        """ Test that the full cross correlation has shape s1 + s2  - 1 """
        im1 = np.random.random(size=(32, 32, 5))
        im2 = np.random.random(size=(32, 32, 5))
        ret = xcorr(im1, im2, mode="full", axes=(0, 1, 2))
        self.assertTupleEqual(ret.shape, (63, 63, 9))

    def test_against_scipy_float(self):
        """ Test against scipy.signal.correlate with float inputs """
        im1 = np.random.random(size=(32, 32, 5))
        im2 = np.random.random(size=(32, 32, 5))

        with self.subTest("mode = 'full'"):
            from_skued = xcorr(im1, im2, mode="full", axes=(0, 1, 2))
            from_scipy = correlate(im1, im2, mode="full")

            self.assertTrue(np.allclose(from_scipy, from_skued))

        with self.subTest("mode = 'same'"):
            from_skued = xcorr(im1, im2, mode="same", axes=(0, 1, 2))
            from_scipy = correlate(im1, im2, mode="same")

            self.assertTrue(np.allclose(from_scipy, from_skued))

    def test_against_scipy_complex(self):
        """ Test against scipy.signal.correlate with complex inputs """
        im1 = np.random.random(size=(32, 32, 5)) + 1j * np.random.random(
            size=(32, 32, 5)
        )
        im2 = np.random.random(size=(32, 32, 5)) + 1j * np.random.random(
            size=(32, 32, 5)
        )

        with self.subTest("mode = 'full'"):
            from_skued = xcorr(im1, im2, mode="full", axes=(0, 1, 2))
            from_scipy = correlate(im1, im2, mode="full")

            self.assertTrue(np.allclose(from_scipy, from_skued))

        with self.subTest("mode = 'same'"):
            from_skued = xcorr(im1, im2, mode="same", axes=(0, 1, 2))
            from_scipy = correlate(im1, im2, mode="same")

            self.assertTrue(np.allclose(from_scipy, from_skued))


class TestMNXC(unittest.TestCase):
    def test_autocorrelation(self):
        """Masked normalized cross-correlation between identical arrays
		should reduce to an autocorrelation even with random masks."""
        # See random number generator for reproducible results
        np.random.seed(23)

        arr1 = camera()
        arr2 = camera()

        # Random masks with 75% of pixels being valid
        m1 = np.random.choice([True, False], arr1.shape, p=[3 / 4, 1 / 4])
        m2 = np.random.choice([True, False], arr2.shape, p=[3 / 4, 1 / 4])

        xcorr = mnxc(arr1, arr1, m1, m1, axes=(0, 1), mode="same", overlap_ratio=0).real
        max_index = np.unravel_index(np.argmax(xcorr), xcorr.shape)

        # Autocorrelation should have maximum in center of array
        self.assertAlmostEqual(xcorr.max(), 1)
        self.assertTrue(np.allclose(max_index, np.array(arr1.shape) / 2))

    def test_over_axes(self):
        """Masked normalized cross-correlation over axes should be
		equivalent to a loop over non-transform axes."""
        # See random number generator for reproducible results
        np.random.seed(23)

        arr1 = np.random.random((8, 8, 5))
        arr2 = np.random.random((8, 8, 5))

        m1 = np.random.choice([True, False], arr1.shape)
        m2 = np.random.choice([True, False], arr2.shape)

        # Loop over last axis
        with_loop = np.empty_like(arr1, dtype=np.complex)
        for index in range(arr1.shape[-1]):
            with_loop[:, :, index] = mnxc(
                arr1[:, :, index],
                arr2[:, :, index],
                m1[:, :, index],
                m2[:, :, index],
                axes=(0, 1),
                mode="same",
            )

        over_axes = mnxc(arr1, arr2, m1, m2, axes=(0, 1), mode="same")

        self.assertTrue(np.allclose(with_loop, over_axes))

    def test_side_effects(self):
        """Masked normalized cross-correlation should not modify the inputs."""
        shape1 = (2, 2, 2)
        shape2 = (2, 2, 2)

        arr1 = np.zeros(shape1)
        arr2 = np.zeros(shape2)

        # Trivial masks
        m1 = np.ones_like(arr1)
        m2 = np.ones_like(arr2)

        for arr in (arr1, arr2, m1, m2):
            arr.setflags(write=False)

            # If arrays are written to, an exception will be raised.
        mnxc(arr1, arr2, m1, m2)

    def test_output_range(self):
        """Masked normalized cross-correlation should return between 1 and -1."""
        # See random number generator for reproducible results
        np.random.seed(23)

        # Array dimensions must match along non-transformation axes, in
        # this case
        # axis 0
        shape1 = (15, 4, 5)
        shape2 = (15, 12, 7)

        # Initial array ranges between -5 and 5
        arr1 = 10 * np.random.random(shape1) - 5
        arr2 = 10 * np.random.random(shape2) - 5

        # random masks
        m1 = np.random.choice([True, False], arr1.shape)
        m2 = np.random.choice([True, False], arr2.shape)

        xcorr = mnxc(arr1, arr2, m1, m2, axes=(1, 2))

        self.assertLessEqual(xcorr.max(), 1)
        self.assertGreaterEqual(xcorr.min(), -1)

    def test_mismatched_dimensions(self):
        """Masked normalized cross-correlation should raise an error if array
		dimensions along non-transformation axes are mismatched."""
        shape1 = (23, 1, 1)
        shape2 = (6, 2, 2)

        arr1 = np.zeros(shape1)
        arr2 = np.zeros(shape2)

        # Trivial masks
        m1 = np.ones_like(arr1)
        m2 = np.ones_like(arr2)

        with self.assertRaises(ValueError):
            mnxc(arr1, arr2, m1, m2, axes=(1, 2))

    def test_output_shape(self):
        """Masked normalized cross-correlation should return a shape
		of N + M + 1 for each transform axis in 'full' mode, and unchanged dimensions 
		in 'same' mode."""
        shape1 = (15, 4, 5)
        shape2 = (6, 12, 7)
        expected_full_shape = tuple(np.array(shape1) + np.array(shape2) - 1)
        expected_same_shape = shape1

        arr1 = np.zeros(shape1)
        arr2 = np.zeros(shape2)
        # Trivial masks
        m1 = np.ones_like(arr1)
        m2 = np.ones_like(arr2)

        full_xcorr = mnxc(arr1, arr2, m1, m2, axes=(0, 1, 2), mode="full")
        self.assertTupleEqual(full_xcorr.shape, expected_full_shape)

        same_xcorr = mnxc(arr1, arr2, m1, m2, axes=(0, 1, 2), mode="same")
        self.assertTupleEqual(same_xcorr.shape, expected_same_shape)


if __name__ == "__main__":
    unittest.main()
