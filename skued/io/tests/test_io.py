# -*- coding: utf-8 -*-
import unittest
from pathlib import Path

import numpy as np

from .. import diffread

TEST_MIB = Path(__file__).parent / 'test.mib'

class TestDiffRead(unittest.TestCase):

    def test_on_merlin_image_binary(self):
        """ Test diffread() on Merlin Image Binary (.mib) """
        im = diffread(TEST_MIB)
        self.assertEqual(im.shape, (256, 256))
        self.assertEqual(im.dtype,  np.dtype('>u2'))

if __name__ == '__main__':
    unittest.main()
