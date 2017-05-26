# -*- coding: utf-8 -*-
from .. import Crystal
import unittest

import os

import numpy as np

BASE_PATH = os.path.join(os.path.abspath('skued'), 'structure', 'tests', 'cifs')
TEST_FILES = list(filter(lambda path: path.endswith('.cif'), os.listdir(BASE_PATH)))

@unittest.skip('')
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