# -*- coding: utf-8 -*-

import unittest
from random import randint, random, seed
from ..thin_films import film_optical_coefficients
from math import sqrt

seed(23)

class TestLightFilmInteraction(unittest.TestCase):

    def test_conservation_power(self):
        """ Test that R + T + A = 1 always """

        wavelength = randint(200, 1000)
        thickness = randint(50, 200)
        n_film = random() + 1j*random()
        n_subs = random() + 1j*random()
        
        R, T, A = film_optical_coefficients(wavelength, thickness, n_film, n_subs)

        self.assertEqual(R + T + A, 1)
    
    def test_correctness(self):
        """ Test for absorption values of VO2 on SiN substrate. Values provided from Martin R. Otto """
        _, _, A = film_optical_coefficients(800, thickness = 90, n_film = 2.9 + 1j*0.43, n_substrate = 1.9962)
        self.assertAlmostEqual(A, 0.312495, places = 5)

if __name__ == '__main__':
    unittest.main()