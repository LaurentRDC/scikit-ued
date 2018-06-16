# -*- coding: utf-8 -*-
import pickle
import unittest
from copy import deepcopy
import collections.abc as abc

import numpy as np
from itertools import islice

from .. import Atom, AtomicStructure, Crystal

class TestAtomicStructure(unittest.TestCase):

    def setUp(self):
        self.substructure = AtomicStructure(atoms = [Atom('U', [0,0,0])])
        self.structure = AtomicStructure(atoms = [Atom('Ag', [0,0,0]), 
                                                  Atom('Ag', [1,1,1])],
                                         substructures = [self.substructure])
    
    def test_iteration(self):
        """ Test iteration of AtomicStructure yields from orphan atoms and substructure atoms alike """
        elements = [atm.element for atm in self.structure]
        self.assertTrue(len(elements), 3)
    
    def test_itersorted(self):
        """ Test that AtomicStructure.itersorted() works as expected """
        sorted_from_structure = list(self.structure.itersorted())
        sorted_from_list      = list(sorted(self.structure, key = lambda a: a.element))

        self.assertListEqual(sorted_from_structure, sorted_from_list)
    
    def test_chemical_composition_trivial(self):
        """ Test that AtomicStructure.chemical_composition works as expected """
        expected = {'U': 1/3, 'Ag':2/3}
        self.assertDictEqual(self.structure.chemical_composition, expected)
        
    def test_chemical_composition_add_to_unity(self):
        """ Test that AtomicStructure.chemical_composition always adds up to 1 """
        # Faster to create a large atomic structure from a Crystal object
        # Testing for 10 crystal structures only
        for name in islice(Crystal.builtins, 10):
            structure = AtomicStructure(atoms = Crystal.from_database(name))
            self.assertAlmostEqual(sum(structure.chemical_composition.values()), 1)

    def test_length(self):
        """ Test the __len__ methods """
        self.assertTrue(len(self.structure), 3)
    
    def test_containership_substructures(self):
        """ Test that containership works on substructure and atoms separately """
        self.assertIn(self.substructure, self.structure)
        self.assertNotIn(self.structure, self.substructure)
    
    def test_containership_atoms(self):
        """ Test that atom containership testing is working, even in substructures """
        atm = next(iter(self.substructure))
        self.assertIn(atm, self.structure)
    
    def test_equality(self):
        """ Test that AtomicStructure is equal to itself but not others """
        self.assertEqual(self.structure, self.structure)
        self.assertEqual(self.structure, deepcopy(self.structure))
        self.assertNotEqual(self.structure, self.substructure)

    def test_array(self):
        """ Test AtomicStructure.__array__ """
        arr = np.array(self.structure)
        self.assertSequenceEqual(arr.shape, (len(self.structure), 4))

    def test_picklable(self):
        """ Test that Crystal instances can be pickled, and that the unpickled instance
        is identical to the source """
        pickled = pickle.dumps(self.structure)
        unpickled = pickle.loads(pickled)
        self.assertEqual(self.structure, unpickled)
    
    def test_abstract_base_classes(self):
        """ Test that AtomicStructure fits with collections.abc module """
        for abstract_base_class in (abc.Hashable, abc.Iterable, abc.Sized):
            self.assertIsInstance(self.structure, abstract_base_class)

if __name__ == '__main__':
    unittest.main()