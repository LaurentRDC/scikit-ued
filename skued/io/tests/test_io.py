# -*- coding: utf-8 -*-
import os
import unittest
from datetime import datetime
from pathlib import Path

import numpy as np
from skimage.io import imsave

from .. import diffread, dmread, imibread, mibheader, mibread
from ...utils import suppress_warnings

TEST_MIB = Path(__file__).parent / 'test.mib'
TEST_MIB_MULTI = Path(__file__).parent / 'test_multi.mib'

TEST_DM3 = Path(__file__).parent / 'bismuth.dm3'
TEST_DM4 = Path(__file__).parent / 'bismuth.dm4'

TEST_PNG = Path(__file__).parent / 'png_test.png'

class TestDiffRead(unittest.TestCase):

    def test_on_merlin_image_binary(self):
        """ Test diffread() on Merlin Image Binary (.mib) """
        im = diffread(TEST_MIB)
        self.assertEqual(im.shape, (256, 256))
        self.assertEqual(im.dtype,  np.dtype('>u2'))
    
    def test_on_dm3_vs_dm4_image(self):
        """ Test that diffread() works on DM3 images """
        im3 = diffread(TEST_DM3)
        im4 = diffread(TEST_DM4)
        
        with self.subTest('DM3'):
            self.assertEqual(im3.shape, (2048, 2048))
            self.assertEqual(im3.dtype, np.dtype('int8'))

        with self.subTest('DM4'):
            self.assertEqual(im4.shape, (2048, 2048))
            self.assertEqual(im4.dtype, np.dtype('int8'))

        with self.subTest('DM3 vs. DM4'):
            self.assertTrue(np.allclose(im3, im4))

    def test_on_tiff(self):
        """ Test diffread() on tiff files """
        im = np.random.randint(0, 127, size = (512, 512))
        path = Path('.\\test_tif.tif')

        # Annoying low contrast warning
        with suppress_warnings():
            imsave(str(path), im)

        from_skued = diffread(path)
        self.assertTrue(np.allclose(im, from_skued))
        os.remove(path)
    
    def test_on_skimage_png(self):
        """ Test the last resort of using skimage.io for pngs """
        from_skimage = diffread(TEST_PNG)

        self.assertTupleEqual(from_skimage.shape, (256,256))
        self.assertTrue(np.allclose(from_skimage, np.ones_like(from_skimage)))

class TestMIBHeader(unittest.TestCase):

    def test_header(self):
        """ Test that header parsing of MIB files is working as intended """
        header = mibheader(TEST_MIB)
        
        true_value = {'ID'       :'MQ1',
                      'seq_num'  : 1,
                      'offset'   : 384,
                      'nchips'   : 1,
                      'shape'   : ( 256, 256 ),
                      'dtype'    : np.dtype('>u2'),
                      'timestamp': datetime(2018, 1, 19, 20, 55, 10, 966026).timestamp()}

        self.assertDictEqual(header, true_value)

class TestMIBRead(unittest.TestCase):

    
    def test_imibread(self):
        """ Test the generator version of mibread() """
        gen = imibread(TEST_MIB)
        arr = next(gen)
        self.assertEqual(arr.shape, (256, 256))
        self.assertEqual(arr.dtype,  np.dtype('>u2'))

    def test_mibread(self):
        """ Test that the array extracted from a test MIB files has the
        expected attributes """
        arr = mibread(TEST_MIB)
        self.assertEqual(arr.shape, (256, 256))
        self.assertEqual(arr.dtype, np.dtype('>u2'))
    
    def test_mibread_multi(self):
        """ Test that the array extracted from a test MIB file containing
        multiple images has the expected attributes """
        arr = mibread(TEST_MIB_MULTI)
        self.assertEqual(arr.shape, (256, 256, 9))
        self.assertEqual(arr.dtype, np.dtype('>u1'))

if __name__ == '__main__':
    unittest.main()
