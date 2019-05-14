# -*- coding: utf-8 -*-
import unittest
import numpy as np
from random import randint

from skued import structure_factor, affe
from crystals import Crystal, Atom


class TestElectronFormFactor(unittest.TestCase):
    def test_side_effects(self):
        nG = np.random.random(size=(16, 32))
        nG.setflags(write=False)  # if nG is written to, Exception is raised
        affe(Atom("He", coords=[0, 0, 0]), nG)

    def test_out_shape(self):
        nG = np.random.random(size=(16, 32))
        eff = affe(Atom("He", coords=[0, 0, 0]), nG)
        self.assertSequenceEqual(eff.shape, nG.shape)

    def test_int(self):
        """ Test that affe(int, ...) also works """
        atomic_number = randint(1, 103)
        nG = np.random.random(size=(16, 32))

        from_int = affe(atomic_number, nG)
        from_atom = affe(Atom(atomic_number, [0, 0, 0]), nG)

        self.assertTrue(np.allclose(from_int, from_atom))

    def test_str(self):
        """ Test that affe(str, ...) also works """
        # Try with Chlorine (Z = 17)
        atomic_number = 17
        nG = np.random.random(size=(16, 32))

        from_int = affe(atomic_number, nG)
        from_str = affe("Cl", nG)

        self.assertTrue(np.allclose(from_int, from_str))


class TestStructureFactor(unittest.TestCase):
    def setUp(self):
        self.crystal = Crystal.from_database(next(iter(Crystal.builtins)))

    def test_shape_and_dtype(self):
        """ Test that output of structure_factor is same shape as input,
        and that the dtype is complex """
        h, k, l = np.meshgrid([1, 2, 3], [1, 2, 3], [1, 2, 3])
        sf = structure_factor(self.crystal, h, k, l)

        self.assertSequenceEqual(sf.shape, h.shape)
        self.assertEqual(sf.dtype, np.complex)


if __name__ == "__main__":
    unittest.main()
