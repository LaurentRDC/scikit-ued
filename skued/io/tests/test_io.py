# -*- coding: utf-8 -*-
import unittest
from pathlib import Path
from skimage.io import imsave
import os
import numpy as np

from .. import diffread

TEST_MIB = Path(__file__).parent / 'test.mib'

class TestDiffRead(unittest.TestCase):

    def test_on_merlin_image_binary(self):
        """ Test diffread() on Merlin Image Binary (.mib) """
        im = diffread(TEST_MIB)
        self.assertEqual(im.shape, (256, 256))
        self.assertEqual(im.dtype,  np.dtype('>u2'))

    def test_on_tiff(self):
        """ Test diffread() on tiff files """
        im = np.random.randint(0, 127, size = (512, 512))
        path = Path('.\\test_tif.tif')
        imsave(str(path), im)

        from_skued = diffread(path)
        self.assertTrue(np.allclose(im, from_skued))
        os.remove(path)

if __name__ == '__main__':
    unittest.main()
