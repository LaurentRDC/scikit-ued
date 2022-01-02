# -*- coding: utf-8 -*-


from random import randint, random, seed

from skued import film_optical_coefficients

seed(23)


def test_film_interaction_conservation_power():
    """Test that R + T + A = 1 always"""

    wavelength = randint(200, 1000)
    thickness = randint(50, 200)
    n_film = random() + 1j * random()
    n_subs = random() + 1j * random()

    R, T, A = film_optical_coefficients(wavelength, thickness, n_film, n_subs)

    assert R + T + A == 1


def test_film_interaction_correctness():
    """Test for absorption values of VO2 on SiN substrate. Values provided from Martin R. Otto"""
    _, _, A = film_optical_coefficients(
        800, thickness=90, n_film=2.9 + 1j * 0.43, n_substrate=1.9962
    )
    assert round(abs(A - 0.312495), 5) == 0
