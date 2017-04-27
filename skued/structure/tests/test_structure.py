
from copy import deepcopy as copy
import numpy as np
from numpy.linalg import norm
from .. import Atom, Structure, AtomicStructure, CompositeStructure

import unittest

np.random.seed(23)

class TestStructure(unittest.TestCase):

    def setUp(self):
        self.structure = Structure(items = [Atom('C', coords = [1,0,0]), 
                                            Atom('Cu', coords = [0,1,0]),
                                            Atom('O', coords = [0,0,1])])
    
    def test_sort(self):
        """ Test sorting by atomic weight """
        key = lambda atom: atom.weight

        items = copy(self.structure.items)
        items.sort(key = key)

        self.structure.sort(attribute = 'weight')
        self.assertSequenceEqual([item.element for item in items], 
                                 [s.element for s in self.structure.items])
    
    def test_find_item(self):
        """ Test finding an item within structure items """
        copper = self.structure.find_item('Cu', key = 'element')
        self.assertEqual(copper.element, 'Cu')

    def test_transform(self):
        pass
    
    def test_rotate(self):
        pass
    
    def translate(self):
        pass