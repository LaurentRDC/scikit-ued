# -*- coding: utf-8 -*-
import unittest
import numpy as np
from numpy.linalg import norm
from .. import Lattice, lattice_vectors_from_parameters

np.random.seed(23)

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
        self.assertAlmostEqual(norm(a1), 1)
        self.assertAlmostEqual(norm(a2), 2)
        self.assertAlmostEqual(norm(a3), 3)

    def test_triclinic(self):
        """ alpha, beta, gama =\= 90 """
        a1, a2, a3 = lattice_vectors_from_parameters(1, 2, 3, 75, 40, 81)
        self.assertAlmostEqual(norm(a1), 1)
        self.assertAlmostEqual(norm(a2), 2)
        self.assertAlmostEqual(norm(a3), 3)
    
    def test_reciprocal_and_back(self):
        """ Create lattice from parameters, take reciprocal twice, 
        and see if the parameters have changed. """
        triclinic = (3, 4, 20, 45, 90, 126)
        triclinic2 = Lattice.from_parameters(*triclinic).reciprocal.reciprocal.lattice_parameters
        self.assertTrue(np.allclose(triclinic, triclinic2))

if __name__ == '__main__':
    unittest.main()