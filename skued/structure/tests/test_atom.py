# -*- coding: utf-8 -*-
from random import choice, seed
import numpy as np
from .. import Atom, atomic_number
import unittest
from copy import deepcopy
seed(23)

ELEMENTS = list(atomic_number.keys())

class TestAtom(unittest.TestCase):
    def setUp(self):
        self.atom = a = Atom(choice([e for e in ELEMENTS if atomic_number[e] < 104]), 
                             np.random.random((3,)))
    
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