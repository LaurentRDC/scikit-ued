# -*- coding: utf-8 -*-
from math import radians, degrees
import unittest
import numpy as np
from numpy.linalg import norm
from ... import rotation_matrix
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
        self.lattice = Lattice(np.eye(3))
    
    def test_euclidian_lattice(self):
        self.assertTrue(np.array_equal(self.lattice.a1, [1,0,0]))
        self.assertTrue(np.array_equal(self.lattice.a2, [0,1,0]))
        self.assertTrue(np.array_equal(self.lattice.a3, [0,0,1]))
    
    def test_volume(self):
        self.assertEqual(self.lattice.volume, 1)

class TestLatticeTransform(unittest.TestCase):
    
    def setUp(self):
        self.lattice = Lattice(np.random.random((3,3)))

    def test_trivial_transformation(self):
        before = np.array(self.lattice.lattice_vectors, copy = True)
        self.lattice.transform(np.eye(3))
        after = np.array(self.lattice.lattice_vectors, copy = True)

        self.assertSequenceEqual(tuple(before.ravel()), tuple(after.ravel()))
    
    def test_wraparound_rotation(self):
        """ Test that a rotation by 360 degrees yields the same lattice """
        before = np.array(self.lattice.lattice_vectors, copy = True)
        self.lattice.transform(rotation_matrix(radians(360), axis = np.random.random((3,))))
        after = np.array(self.lattice.lattice_vectors, copy = True)

        for x1, x2 in zip(tuple(before.ravel()), tuple(after.ravel())):
            self.assertAlmostEqual(x1, x2)
    
    def test_transform_back_and_forth(self):
        before = np.array(self.lattice.lattice_vectors, copy = True)
        transf = np.random.random((3,3))
        self.lattice.transform(transf)
        self.lattice.transform(np.linalg.inv(transf))
        after = np.array(self.lattice.lattice_vectors, copy = True)

        for x1, x2 in zip(tuple(before.ravel()), tuple(after.ravel())):
            self.assertAlmostEqual(x1, x2)

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

class TestLatticeMillerScattering(unittest.TestCase):

    def setUp(self):
        self.lattice = Lattice(np.random.random((3,3)))

    def test_scattering_vector_trivial(self):
        """ Test that Lattice.scattering_vectors is working """
        self.assertSequenceEqual(self.lattice.scattering_vector(0,0,0), (0,0,0))
    
    def test_miller_indices_trivial(self):
        """ Test that Lattice.miller_indices is working """
        self.assertSequenceEqual(self.lattice.miller_indices(0,0,0), (0,0,0))
    
    def test_back_and_forth(self):
        """ Test that Lattice.miller_indices and Lattice.scattering_vector are
        reciprocal to each other """
        h, k, l = np.random.randint(-10, 10, size = (3,))
        Gx, Gy, Gz = self.lattice.scattering_vector(h, k, l)
        hp, kp, lp = self.lattice.miller_indices(Gx, Gy, Gz)
        
        self.assertAlmostEqual(h, float(hp))
        self.assertAlmostEqual(k, float(kp))
        self.assertAlmostEqual(l, float(lp))

if __name__ == '__main__':
    unittest.main()