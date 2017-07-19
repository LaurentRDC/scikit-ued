# -*- coding: utf-8 -*-
from math import radians
from copy import deepcopy, copy
from random import choice, seed
from itertools import permutations
import numpy as np
from .. import Crystal, Atom, Lattice, graphite
from ... import rotation_matrix, transform
import unittest

#seed(23)

class TestBoundedReflections(unittest.TestCase):

    def setUp(self):
        self.crystal = Crystal.from_database(next(iter(Crystal.builtins)))

    def test_bounded_reflections_negative(self):
        """ Test that negative reflection bounds raise an Exception.
        Otherwise, an infinite number of reflections will be generated """
        with self.assertRaises(ValueError):
            hkl = list(self.crystal.bounded_reflections(-1))
    
    def test_bounded_reflections_zero(self):
        """ Check that bounded_reflections returns (000) for a zero bound """
        h, k, l = self.crystal.bounded_reflections(0)
        [self.assertEqual(len(i), 1) for i in (h, k, l)]
        [self.assertEqual(i[0], 0) for i in (h, k, l)]
    
    def test_bounded_reflections_all_within_bounds(self):
        """ Check that every reflection is within the bound """
        bound = 10
        Gx, Gy, Gz = self.crystal.scattering_vector(*self.crystal.bounded_reflections(nG = bound))
        norm_G = np.sqrt(Gx**2 + Gy**2 + Gz**2)
        self.assertTrue(np.all(norm_G <= bound))

class TestCrystalRotations(unittest.TestCase):

    def setUp(self):
        self.crystal = Crystal.from_database(next(iter(Crystal.builtins)))
    
    def test_crystal_equality(self):
        """ Tests that Crystal.__eq__ is working properly """
        self.assertEqual(self.crystal, self.crystal)

        cryst2 = deepcopy(self.crystal)
        cryst2.transform(2*np.eye(3)) # This stretches lattice vectors, symmetry operators
        self.assertFalse(self.crystal is cryst2)
        self.assertNotEqual(self.crystal, cryst2)

        cryst2.transform(0.5*np.eye(3))
        self.assertEqual(self.crystal, cryst2)
    
    def test_trivial_rotation(self):
        """ Test rotation by 360 deg around all axes. """
        unrotated = deepcopy(self.crystal)
        r = rotation_matrix(radians(360), [0,0,1])
        self.crystal.transform(r)

        self.assertEqual(self.crystal, unrotated)
    
    def test_identity_transform(self):
        """ Tests the trivial identity transform """
        transf = deepcopy(self.crystal)
        transf.transform(np.eye(3))
        self.assertEqual(self.crystal, transf)
    
    def test_one_axis_rotation(self):
        """ Tests the crystal orientation after rotations. """
        unrotated = deepcopy(self.crystal)
        self.crystal.transform(rotation_matrix(radians(37), [0,1,0]))
        self.assertNotEqual(unrotated, self.crystal)
        self.crystal.transform(rotation_matrix(radians(-37), [0,1,0]))
        self.assertEqual(unrotated, self.crystal)

    def test_wraparound_rotation(self):
        cryst1 = deepcopy(self.crystal)
        cryst2 = deepcopy(self.crystal)

        cryst1.transform(rotation_matrix(radians(22.3), [0,0,1]))
        cryst2.transform(rotation_matrix(radians(22.3 - 360), [0,0,1]))
        self.assertEqual(cryst1, cryst2)
    
class TestCrystalConstructors(unittest.TestCase):

    def test_builtins(self):
        """ Test that all names in Crystal.builtins build without errors """
        for name in Crystal.builtins:
            c = Crystal.from_database(name)
    
    def test_builtins_wrong_name(self):
        """ Test that a name not in Crystal.builtins will raise a ValueError """
        with self.assertRaises(ValueError):
            c = Crystal.from_database('___')
    
    def test_from_cod(self):
        """ Test building a Crystal object from the COD """
        # revision = None and latest revision should give the same Crystal
        c = Crystal.from_cod(1521124)
        c2 = Crystal.from_cod(1521124, revision = 176429)

        self.assertSetEqual(set(c), set(c2))

if __name__ == '__main__':
    unittest.main()