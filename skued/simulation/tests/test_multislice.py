# -*- coding: utf-8 -*-
from .. import sim_mesh
from ... import Crystal
import numpy as np
import unittest
from random import seed, choice

seed(23)

class TestSimMesh(unittest.TestCase):

    # Only parse the CIFs once
    crystals = [Crystal.from_database(name) for name in iter(Crystal.builtins)]

    def test_sim_mesh_range(self):
        """ Test that the spatial meshes returned by sim_mesh() fit
        an integer amount of unit cells """

        for crystal in self.crystals:
            with self.subTest('Crystal = {}'.format(crystal.source)):
                X, Y, _, _ = sim_mesh(crystal, resolution = (512, 512))
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
                X, Y, _, _ = sim_mesh(crystal, resolution = (512, 512))
                per_x, per_y, _ = crystal.periodicity
                dx, dy = X[1,1] - X[0,0], Y[1,1] - Y[0,0]

                # Rounding down to 5 digits because of floating errors
                ndx = round(per_x / dx, ndigits = 8)
                ndy = round(per_y / dy, ndigits = 8)

                self.assertAlmostEqual(ndx, int(ndx))
                self.assertAlmostEqual(ndy, int(ndy))
                
if __name__ == '__main__':
	unittest.main()