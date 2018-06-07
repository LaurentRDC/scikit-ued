# -*- coding: utf-8 -*-
import unittest

import numpy as np
from .. import Crystal, plane_mesh, potential_map, powdersim
from ..utils import suppress_warnings

class TestPotentialMap(unittest.TestCase):

    def setUp(self):
        self.crystal = Crystal.from_database('vo2-rutile')
        self.q = np.linspace(1, 5, 256)
        self.I = powdersim(self.crystal, self.q)
    
    def test_shape(self):
        """ Test that potential_map returns a map with the same shape as the mesh """
        aR1, aR2, aR3 = self.crystal.lattice_vectors
        extent = np.arange(0, 5, 0.1)

        with suppress_warnings():
            plane = plane_mesh(aR3, aR1 + aR2, x1 = extent)
        
        potmap = potential_map(self.q, self.I, self.crystal, plane)

        xx, yy, zz = plane
        for arr in plane:
            self.assertTupleEqual(potmap.shape, arr.shape)
        
    def test_positive_intensity(self):
        """ Test that potential_map raises an error if diffraction intensity is not positive """
        aR1, aR2, aR3 = self.crystal.lattice_vectors
        extent = np.arange(0, 5, 0.1)
        with suppress_warnings():
            plane = plane_mesh(aR3, aR1 + aR2, x1 = extent)
        
        self.I[0] = -1
        with self.assertRaises(ValueError):
            potmap = potential_map(self.q, self.I, self.crystal, plane)

    def test_trivial(self):
        """ Test that potential_map calculated from zero intensity is zero everywhere """
        aR1, aR2, aR3 = self.crystal.lattice_vectors
        extent = np.arange(0, 10, 0.1)
        
        with suppress_warnings():
            plane = plane_mesh(aR3, aR1 + aR2, x1 = extent)
        
        potmap = potential_map(self.q, np.zeros_like(self.I), self.crystal, plane)

        self.assertTrue(np.allclose(potmap, 0))

if __name__ == '__main__':
    unittest.main()
