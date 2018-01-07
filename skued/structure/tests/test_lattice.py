# -*- coding: utf-8 -*-
import unittest
from copy import deepcopy
from math import degrees, radians

import numpy as np
from numpy.linalg import norm

from .. import Lattice, LatticeSystem, Crystal
from ... import rotation_matrix

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
    
    def test_equality(self):
        """ Test equality between identical Lattice instances and copies """
        self.assertEqual(self.lattice, self.lattice)
        self.assertEqual(self.lattice, deepcopy(self.lattice))
        
        self.assertNotEqual(self.lattice, Lattice(2*np.eye(3)))

class TestLatticeArray(unittest.TestCase):
    
    def test_shape(self):
        """ Test that array(Lattice(...)) is always 3x3 """
        arr = np.random.random(size = (3,3))
        lattice = Lattice(arr)
        self.assertTrue(np.allclose(arr, np.array(lattice)))
    
    def test_dtype(self):
        """ Test that the data-type of array(Lattice(...)) is respected """
        arr = np.random.random(size = (3,3))
        lattice = Lattice(arr)
        self.assertTrue(np.array(lattice, dtype = np.int).dtype, np.int)

class TestLatticeMeshes(unittest.TestCase):

    def test_frac_mesh(self):
        """ Test that Lattice.frac_mesh is working as expected compared to numpy.meshgrid """
        lattice = Lattice(np.eye(3))
        x = np.linspace(0, 1, num = 8)

        for out, n in zip(lattice.frac_mesh(x), np.meshgrid(x,x,x)):
            self.assertTrue(np.allclose(out, n))

    def test_frac_mesh_two_arr(self):
        """ Test that Lattice.frac_mesh is raising an exception for two input vectors """
        lattice = Lattice(np.eye(3))
        x = np.linspace(0, 1, num = 2)

        with self.assertRaises(ValueError):
            lattice.frac_mesh(x, x)
    
    def test_real_mesh_trivial(self):
        """ Test that Lattice.mesh works identically to Lattice.frac_mesh for trivial lattice """
        lattice = Lattice(np.eye(3))
        x = np.linspace(0, 1, num = 8)

        for frac, real in zip(lattice.frac_mesh(x), lattice.mesh(x)):
            self.assertTrue(np.allclose(frac, real))
    
    def test_real_mesh(self):
        """ Test that Lattice.mesh works as expected """
        lattice = Lattice(2*np.eye(3))
        x = np.linspace(0, 1, num = 8)
        
        # since lattice is a stretched euclidian lattice, we expect
        # a maximum length of 2
        self.assertTrue(np.max(lattice.mesh(x)[0]) == 2)

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
        a1, a2, a3 = Lattice.from_parameters(2,1,5,90,90,90).lattice_vectors
        self.assertTrue(np.allclose(a1, [2,0,0]))
        self.assertTrue(np.allclose(a2, [0,1,0]))
        self.assertTrue(np.allclose(a3, [0,0,5]))

    def test_monoclinic(self):
        """ beta =\= 90 """
        a1, a2, a3 = Lattice.from_parameters(1,2,3, 90, 120, 90).lattice_vectors
        
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
        a1, a2, a3 = Lattice.from_parameters(1, 2, 3, 75, 40, 81).lattice_vectors

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

class TestLatticeSystems(unittest.TestCase):

    def test_cubic_lattice_system(self):
        """ Test that the Lattice.lattice_system attribute is working properly for cubic lattice """
        self.assertEqual(Lattice(2*np.eye(3)).lattice_system,  LatticeSystem.cubic)

    def test_tetragonal_lattice_system(self):
        """ Test that the Lattice.lattice_system attribute is working properly for tetragonal lattice """
        parameters = (2,2,3, 90, 90, 90)
        l = Lattice.from_parameters(*parameters)
        self.assertEqual(l.lattice_system, LatticeSystem.tetragonal)

    def test_rhombohedral_lattice_system(self):
        """ Test that the Lattice.lattice_system attribute is working properly for rhombohedral lattice """
        parameters = (1, 1, 1, 87, 87, 87)
        l = Lattice.from_parameters(*parameters)
        self.assertEqual(l.lattice_system, LatticeSystem.rhombohedral)

    def test_monoclinic_lattice_system(self):
        """ Test that the Lattice.lattice_system attribute is working properly for monoclinic lattice
        including all possible permutations. """
        with self.subTest('permutation 1'):
            parameters = (1, 2, 3, 90, 115, 90)
            l = Lattice.from_parameters(*parameters)
            self.assertEqual(l.lattice_system, LatticeSystem.monoclinic)
        
        with self.subTest('permutation 2'):
            parameters = (2, 3, 1, 115, 90, 90)
            l = Lattice.from_parameters(*parameters)
            self.assertEqual(l.lattice_system, LatticeSystem.monoclinic)

        with self.subTest('permutation 3'):
            parameters = (3, 1, 2, 90, 90, 115)
            l = Lattice.from_parameters(*parameters)
            self.assertEqual(l.lattice_system, LatticeSystem.monoclinic)

    def test_hexagonal_lattice_system(self):
        """ Test that the Lattice.lattice_system attribute is working properly for hexagonal lattice,
        including all possible permutations of lattice parameters. """
        with self.subTest('Gamma = 120deg'):
            parameters = (2,2,3, 90, 90, 120)
            l = Lattice.from_parameters(*parameters)
            self.assertEqual(l.lattice_system, LatticeSystem.hexagonal)
        
        with self.subTest('alpha = 120deg'):
            parameters = (3,2,2, 120, 90, 90)
            l = Lattice.from_parameters(*parameters)
            self.assertEqual(l.lattice_system, LatticeSystem.hexagonal)

        with self.subTest('beta = 120deg'):
            parameters = (2,3,2, 90, 120, 90)
            l = Lattice.from_parameters(*parameters)
            self.assertEqual(l.lattice_system, LatticeSystem.hexagonal)
        
        with self.subTest('Equal lengths'):
            parameters = (2,2,2, 90, 120, 90)
            l = Lattice.from_parameters(*parameters)
            self.assertEqual(l.lattice_system, LatticeSystem.hexagonal)
    
    def test_triclinic_lattice_system(self):
        """ Test that the Lattice.lattice_system attribute is working properly for triclinic lattice """
        l = Lattice.from_parameters(1, 2, 3, 75, 40, 81)
        self.assertEqual(l.lattice_system, LatticeSystem.triclinic)

    def test_graphite(self):
        """ Test that the builtin Crystal for graphite has a hexagonal lattice system """
        graphite = Crystal.from_database('C')
        self.assertEqual(graphite.lattice_system, LatticeSystem.hexagonal)
    
    def test_lead(self):
        """ Test that the builtin Crystal for lead has a cubic lattice system """
        pb = Crystal.from_database('Pb')
        self.assertEqual(pb.lattice_system,  LatticeSystem.cubic)
    
    def test_vo2(self):
        """ Test that the builtin Crystal for monoclinic M1 VO2 has a monoclinic lattice system """
        vo2 = Crystal.from_database('vo2-m1')
        self.assertEqual(vo2.lattice_system, LatticeSystem.monoclinic)


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
