# -*- coding: utf-8 -*-
from random import choice, seed, random, randint
import numpy as np
from .. import Atom, Lattice
from ... import rotation_matrix
import unittest
from copy import deepcopy

seed(23)
np.random.seed(23)

try:
    import ase
    ASE = True
except ImportError:
    ASE = False

def random_transform():
    return rotation_matrix(random(), axis = np.random.random((3,)))

@unittest.skipIf(not ASE, 'ASE not installed or importable')
class TestASEAtom(unittest.TestCase):

    def setUp(self):
        self.atom =  Atom(randint(1, 103), coords = np.random.random((3,)))
    
    def test_back_and_forth(self):
        """ Test that conversion from skued.Atom and ase.Atom is working """
        to_ase = self.atom.ase_atom()
        atom2 = Atom.from_ase(to_ase)
        self.assertEqual(self.atom, atom2)

class TestAtom(unittest.TestCase):

    def setUp(self):
        self.atom =  Atom(randint(1, 103), coords = np.random.random((3,)))
    
    def test_init(self):
        """ Test that Atom can be instantiated with an element str or atomic number """
        by_element = Atom('C', coords = (0,0,0))
        by_number = Atom(6, coords = (0,0,0))
        self.assertEqual(by_element, by_number)
    
    def test_equality(self):
        """ Test __eq__ for atoms """
        other = deepcopy(self.atom)
        self.assertEqual(self.atom, self.atom)
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

if __name__ == '__main__':
    unittest.main()