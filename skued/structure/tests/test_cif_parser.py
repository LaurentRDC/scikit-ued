# -*- coding: utf-8 -*-
from .. import Crystal
import unittest

from os.path import dirname, join

import numpy as np

TEST_FILE = join('skued', 'structure', 'tests', 'cu.cif')

class TestCifParser(unittest.TestCase):

    def setUp(self):
        self.crystal = Crystal.from_cif(TEST_FILE)
    
    def test_atoms(self):
        """ Test that the only atoms are Cu"""
        self.assertIn('Cu', map(lambda atm : atm.element, self.crystal))
        self.assertEqual(len(self.crystal.atoms), 1)
    
    def test_lattice_vectors(self):
        """ Test the correctness of CIF lattice vectors """
        self.assertSequenceEqual(list(map(np.linalg.norm, self.crystal.lattice_vectors)), 
                                 [3.59127,3.59127,3.59127])
    
if __name__ == '__main__':
    unittest.main()