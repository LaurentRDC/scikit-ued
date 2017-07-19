# -*- coding: utf-8 -*-
from random import choice, seed, random
import numpy as np
from .. import Atom, atomic_number, Lattice
from ... import rotation_matrix
import unittest
from copy import deepcopy

seed(23)
np.random.seed(23)

ELEMENTS = list(atomic_number.keys())

def random_transform():
    return rotation_matrix(random(), axis = np.random.random((3,)))

class TestAtom(unittest.TestCase):

    def setUp(self):
        # Atoms with Z larger than 104 have no scattering parameters available
        self.atom =  Atom(choice([e for e in ELEMENTS if atomic_number[e] < 104]), 
                          coords = np.random.random((3,)))
    
    def test_electron_form_factor_side_effects(self):
        """ Test that arrays passed to Atom.electron_form_factors are unchanged, which
        has been a problem in the past. """
        nG = np.random.random(size = (256,))
        copied = np.copy(nG)
        _ = self.atom.electron_form_factor(nG)
        self.assertTrue(np.allclose(nG, copied))
    
    def test_equality(self):
        """ Test __eq__ for atoms """
        other = deepcopy(self.atom)
        self.assertEqual(self.atom, other)

        other.coords = self.atom.coords + 1
        self.assertNotEqual(self.atom, other)

    def test_trivial_transform(self):
        """ Test Atom.transform() with the identity """
        before = np.array(self.atom.coords, copy = True)
        self.atom.transform(np.eye(3))
        after = np.array(self.atom.coords, copy = True)

        self.assertSequenceEqual(tuple(before), tuple(after))
    
    def test_transform_back_and_forth(self):
        """ Test Atom.transform() with a random transformation back and forth """
        before = np.array(self.atom.coords, copy = True)

        transf = random_transform()
        self.atom.transform(transf)
        self.atom.transform(np.linalg.inv(transf))
        after = np.array(self.atom.coords, copy = True)

        # No assert sequence almost equal
        for x1, x2 in zip(tuple(before), tuple(after)):
            self.assertAlmostEqual(x1, x2)


