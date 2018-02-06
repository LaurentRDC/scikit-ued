# -*- coding: utf-8 -*-
import unittest
from os.path import dirname, join
import numpy as np

from .. import diffread

# If Py3.6, test compatibility with pathlib
try:
    from pathlib import Path
except ImportError:
    WITH_PATHLIB = False
else:
    WITH_PATHLIB = True

TEST_MIB = join(dirname(__file__), 'test.mib')

class TestDiffRead(unittest.TestCase):

    def test_on_merlin_image_binary(self):
        """ Test diffread() on Merlin Image Binary (.mib) """
        im = diffread(TEST_MIB)
        self.assertEqual(im.shape, (256, 256))
        self.assertEqual(im.dtype,  np.dtype('>u2'))
    
    @unittest.skipIf(not WITH_PATHLIB, 'pathlib not importable (possibly Python version < 3.6)')
    def test_compat_with_pathlib(self):
        """ Test diffread() on pathlib.Path instead of strings """
        im = diffread(Path(TEST_MIB))
        self.assertEqual(im.shape, (256, 256))
        self.assertEqual(im.dtype,  np.dtype('>u2'))

if __name__ == '__main__':
    unittest.main()