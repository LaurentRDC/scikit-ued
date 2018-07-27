# -*- coding: utf-8 -*-
import unittest

import numpy as np

from ..fourier import fft2freq, limit_bandwidth

class TestFFT2Freq(unittest.TestCase):

    def test_shape(self):
        """ Test that the output from fft2freq has the same shape as the input. """
        extent_x = np.linspace(0, 1)
        extent_y = np.linspace(3, 4, num = 47)

        for indexing in {'ij', 'xy'}:
            with self.subTest('Indexing {}'.format(indexing)):
                x, y = np.meshgrid(extent_x, extent_y, indexing = indexing)
                kx, ky = fft2freq(x, y, indexing = indexing)
                self.assertTupleEqual(x.shape, kx.shape)
                self.assertTupleEqual(y.shape, ky.shape)

    def test_indexing_error(self):
        """ Test that fft2freq correctly raises ValueError for invalid indexing. """
        extent_x = np.linspace(0, 1)
        extent_y = np.linspace(3, 4, num = 47)

        with self.assertRaises(ValueError):
                fft2freq([1,2], [1,2], indexing = 'ab')
    
    def test_vs_fftfreq(self):
        """ Test that the results make sense with respect to 1D case """
        extent_x = np.arange(0, 10, step = 0.1)
        extent_y = np.arange(3, 4, step = 0.1)

        for indexing in {'ij', 'xy'}:
            with self.subTest('Indexing {}'.format(indexing)):
                x, y = np.meshgrid(extent_x, extent_y, indexing = indexing)
                kx, ky = fft2freq(x, y, indexing = indexing)

                # same array creates from the 1d case
                kx_1d = np.fft.fftfreq(len(extent_x), d = 0.1)
                ky_1d = np.fft.fftfreq(len(extent_y), d = 0.1)
                kx2, ky2 = np.meshgrid(kx_1d, ky_1d, indexing = indexing)

                self.assertTrue(np.allclose(kx, kx2))
                self.assertTrue(np.allclose(ky, ky2))

class TestLimitBandwidth(unittest.TestCase):

    def test_idempotence(self):
        """ Test that applying limit_bandwidth more than once has no effect """
        im = np.ones((32,32), dtype = np.float)
        kx, ky = np.meshgrid(np.arange(-16, 16), np.arange(-16, 16))
        k = np.hypot(kx, ky)

        limited1 = limit_bandwidth(im, k, k.max()/2)
        limited2 = limit_bandwidth(limited1, k, k.max()/2)

        self.assertTrue(np.allclose(limited1, limited2))


if __name__ == '__main__':
	unittest.main()
