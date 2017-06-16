# -*- coding: utf-8 -*-
import os
import unittest
from warnings import filterwarnings

import numpy as np

from .. import Crystal, graphite
from ... import transform
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

    def test_correctness(self):
        """ Test the correctness of sym_ops """
        with self.subTest('No translation 1'):
            s = sym_ops("+y, +x, +z")
            m = np.array([0,1,0,0,
                          1,0,0,0,
                          0,0,1,0,
                          0,0,0,1]).reshape((4,4))
            self.assertTrue(np.allclose(s, m))
        
        with self.subTest('No translation 2'):
            s = sym_ops("+x-y, +x+y, +z")
            m = np.array([1,-1,0,0,
                          1,1,0,0,
                          0,0,1,0,
                          0,0,0,1]).reshape((4,4))
            self.assertTrue(np.allclose(s, m))


class TestCIRParser(unittest.TestCase):
    """ Test the CIFParser on all CIF files stored herein """

    def _cif_files(self):
        """ Yields cif files included in scikit-ued """
        for root, _, files in os.walk(os.path.join('skued', 'structure')):
            for name in filter(lambda path: path.endswith('.cif'), files):
                yield os.path.join(root, name)

    def test_compatibility(self):
        """ Test the CIFParser on all CIF files stored herein to check build errors"""
        for name in self._cif_files():
            with self.subTest(name.split('\\')[-1]):
                c = Crystal.from_cif(name)
    
    def test_symmetry_operators(self):
        """ Test that the non-translation part of the symmetry_operators is an invertible
        matrix of determinant 1 | -1 """
        for name in self._cif_files():
            with self.subTest(name.split('\\')[-1]):
                with CIFParser(name) as p:
                    for sym_op in p.symmetry_operators():
                        t = sym_op[:3,:3]
                        self.assertAlmostEqual(abs(np.linalg.det(t)), 1)
    
    def test_graphite(self):
        """ Test CIFParser on C.cif and compare to built-in graphite """
        C_path = os.path.join('skued', 'structure', 'cifs', 'periodic_table', 'C.cif')
        c = Crystal.from_cif(C_path)
        
        self.assertEqual(len(c), len(graphite))
        self.assertAlmostEqual(c.volume, graphite.volume, places = -1)

if __name__ == '__main__':
    unittest.main()
