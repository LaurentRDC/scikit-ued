# -*- coding: utf-8 -*-
import os
import unittest
from warnings import filterwarnings

import numpy as np

from .. import Crystal
from ..cif_parser import CIFParser, sym_ops

filterwarnings('ignore', category = UserWarning)

class TestSymOpsParsing(unittest.TestCase):
    """ Test the sym_ops parsing function """

    def test_trivial(self):
        """ Test the identity transform "+x, +y, +z" """
        with self.subTest('string'):
            s = sym_ops("+x, +y, +z")
            self.assertTrue(np.allclose(s, np.eye(4)))
        
        with self.subTest('iterable of strings'):
            s = sym_ops(['+x', '+y', '+z'])
            self.assertTrue(np.allclose(s, np.eye(4)))

class TestParserCompatibility(unittest.TestCase):
    """ Test the CIFParser on all CIF files stored herein """

    def test_all_files(self):
        for root, _, files in os.walk(os.path.join('skued', 'structure', 'tests')):
            for name in filter(lambda path: path.endswith('.cif'), files):
                with self.subTest(name.split('\\')[-1]):
                    full_path = os.path.join(root, name)
                    c = Crystal.from_cif(full_path)

if __name__ == '__main__':
    unittest.main()
