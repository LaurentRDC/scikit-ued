# -*- coding: utf-8 -*-
from .. import sim_mesh, weak_phase, multislice
from ... import Crystal
from itertools import islice
import numpy as np
import unittest
from random import seed

seed(23)

class TestSimMesh(unittest.TestCase):

    # Only parse the CIFs once
    crystals = [Crystal.from_database(name) for name in iter(Crystal.builtins)]

    def test_sim_mesh_range(self):
        """ Test that the spatial meshes returned by sim_mesh() fit
        an integer amount of unit cells """

        for crystal in self.crystals:
            with self.subTest('Crystal = {}'.format(crystal.source)):
                X, Y, _, _ = sim_mesh(crystal, resolution = (256, 256))
                per_x, per_y, _ = crystal.periodicity

                # Rounding down to 5 digits because of floating errors
                ncellx = round((X.max() - X.min()) / per_x, ndigits = 8)
                ncelly = round((Y.max() - Y.min()) / per_y, ndigits = 8)

                self.assertAlmostEqual(ncellx, int(ncellx))
                self.assertAlmostEqual(ncelly, int(ncelly))
    
    def test_sim_mesh_spacing(self):
        """ Test that the spatial spacing of sim_mesh divides a crystal unit cell """
        for crystal in self.crystals:
            with self.subTest('Crystal = {}'.format(crystal.source)):
                X, Y, _, _ = sim_mesh(crystal, resolution = (256, 256))
                per_x, per_y, _ = crystal.periodicity
                dx, dy = X[1,1] - X[0,0], Y[1,1] - Y[0,0]

                # Rounding down to 5 digits because of floating errors
                ndx = round(per_x / dx, ndigits = 8)
                ndy = round(per_y / dy, ndigits = 8)

                self.assertAlmostEqual(ndx, int(ndx))
                self.assertAlmostEqual(ndy, int(ndy))

def wavef_norm(wavefunction):
    return np.sum(np.abs(wavefunction)**2)

class TestWeakPhase(unittest.TestCase):

    # Only test a few Crystal
    crystals = list(islice( (Crystal.from_database(name) for name in iter(Crystal.builtins)), 5))

    def test_intensity_conservation(self):
        """ Test that wavefunction intensity is conserved """
        for crystal in self.crystals:
            with self.subTest('Crystal = {}'.format(crystal.source)):
                X, Y, _, _ = sim_mesh(crystal, resolution = (128, 128))
                initial = np.ones_like(X)
                exit_wave = weak_phase(crystal, 90, initial, resolution = (128, 128))

                self.assertEqual(wavef_norm(initial), wavef_norm(exit_wave))

@unittest.skip('Not ready for primetime')
class TestMultislice(unittest.TestCase):

    # Only test a few Crystal
    crystals = list(islice( (Crystal.from_database(name) for name in iter(Crystal.builtins)), 5))

    def test_output_shape(self):
        """ Test that wavefunction intensity is conserved """
        for crystal in self.crystals:
            with self.subTest('Crystal = {}'.format(crystal.source)):
                X, Y, _, _ = sim_mesh(crystal, resolution = (128, 128))
                initial = np.ones_like(X)
                exit_wave = multislice(crystal, 90, thickness = 10, resolution = (128, 128))

                self.assertSequenceEqual(exit_wave.shape, initial.shape)
if __name__ == '__main__':
	unittest.main()