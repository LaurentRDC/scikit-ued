# -*- coding: utf-8 -*-
from math import degrees
import unittest
import numpy as np
from numpy.linalg import norm
from .. import Lattice, lattice_vectors_from_parameters

np.random.seed(23)

# Angle between vectors 
# https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249
def angle_between(v1, v2):
    """ Returns the angle in degrees between vectors `v1` and `v2`"""
    v1 /= norm(v1)
    v2 /= norm(v2)
    rad = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
    return degrees(rad)

class TestEuclidianLattice(unittest.TestCase):
    
    def setUp(self):
        self.lattice = Lattice()
    
    def test_euclidian_lattice(self):
        self.assertTrue(np.array_equal(self.lattice.a1, [1,0,0]))
        self.assertTrue(np.array_equal(self.lattice.a2, [0,1,0]))
        self.assertTrue(np.array_equal(self.lattice.a3, [0,0,1]))
    
    def test_volume(self):
        self.assertEqual(self.lattice.volume, 1)

class TestLatticeParameters(unittest.TestCase):

    def test_orthorombic(self):
        """ alpha = beta = gamma = 90"""
        a1, a2, a3 = lattice_vectors_from_parameters(2,1,5,90,90,90)
        self.assertTrue(np.allclose(a1, [2,0,0]))
        self.assertTrue(np.allclose(a2, [0,1,0]))
        self.assertTrue(np.allclose(a3, [0,0,5]))

    def test_monoclinic(self):
        """ beta =\= 90 """
        a1, a2, a3 = lattice_vectors_from_parameters(1,2,3, 90, 120, 90)
        
        with self.subTest('Lengths'):
            self.assertAlmostEqual(norm(a1), 1)
            self.assertAlmostEqual(norm(a2), 2)
            self.assertAlmostEqual(norm(a3), 3)
        
        with self.subTest('Angles'):
            self.assertAlmostEqual(angle_between(a2, a3), 90)
            self.assertAlmostEqual(angle_between(a3, a1), 120)
            self.assertAlmostEqual(angle_between(a1, a2), 90)

    def test_triclinic(self):
        """ alpha, beta, gama =\= 90 """
        a1, a2, a3 = lattice_vectors_from_parameters(1, 2, 3, 75, 40, 81)

        with self.subTest('Lengths'):
            self.assertAlmostEqual(norm(a1), 1)
            self.assertAlmostEqual(norm(a2), 2)
            self.assertAlmostEqual(norm(a3), 3)

        with self.subTest('Angles'):
            self.assertAlmostEqual(angle_between(a2, a3), 75)
            self.assertAlmostEqual(angle_between(a3, a1), 40)
            self.assertAlmostEqual(angle_between(a1, a2), 81)
    
    def test_reciprocal_and_back(self):
        """ Create lattice from parameters, take reciprocal twice, 
        and see if the parameters have changed. """
        triclinic = (3, 4, 20, 45, 90, 126)
        triclinic2 = Lattice.from_parameters(*triclinic).reciprocal.reciprocal.lattice_parameters
        self.assertTrue(np.allclose(triclinic, triclinic2))

if __name__ == '__main__':
    unittest.main()