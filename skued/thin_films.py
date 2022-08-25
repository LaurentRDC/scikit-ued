# -*- coding: utf-8 -*-
"""
Calculations on thin films 
==========================
"""
from cmath import cos, exp, pi, sin


def film_optical_coefficients(wavelength, thickness, n_film, n_substrate=0):
    """
    Calculate the reflection, transmission, and absorption coefficients
    of a thin-film (possibly on a substrate).

    Parameters
    ----------
    wavelength : float
        Wavelength of the incident radiation [nm].
    thickness : float
        Thickness of the film [nm].
    n_film : complex or float
        Complex refractive index of the film material.
    n_substrate : complex or float, optional
        Complex refractive index of the substrate material.

    Returns
    -------
    R, T, A : float
        Reflection, transmission, and absorption coefficients, respectively.

    References
    ----------
    .. [#] Tomlin, Brit. J. Appl. Phys. (J. Phys. D) ser. 2. vol. 1 1968
    """
    # Lengths in meters
    wavelength *= 10e-9
    thickness *= 10e-9

    n_film, n_subtrate = complex(n_film), complex(n_substrate)

    # Separate refractive index from absorption
    n0 = 1  # Vacuum
    n1, k1 = n_film.real, n_film.imag
    n2, k2 = n_substrate.real, n_substrate.imag

    # The following are simplifications based on a notebook by Martin R. Otto
    # See also Tomlin, Brit. J. Appl. Phys. (J. Phys. D) ser. 2. vol. 1 1968
    g1 = (n0**2 - n1**2 - k1**2) / ((n0 + n1) ** 2 + k1**2)
    g2 = (n1**2 - n2**2 + k1**2 - k2**2) / ((n1 + n2) ** 2 + (k1 + k2) ** 2)

    h1 = 2 * n0 * k1 / ((n0 + n1) ** 2 + k1**2)
    h2 = 2 * (n1 * k2 - n2 * k1) / ((n1 + n2) ** 2 + (k1 + k2) ** 2)

    alpha1 = 2 * pi * k1 * thickness / wavelength
    gamma1 = 2 * pi * n1 * thickness / wavelength

    A = 2 * (g1 * g2 + h1 * h2)
    B = 2 * (g1 * h2 - g2 * h1)
    C1 = 2 * (g1 * g2 - h1 * h2)
    D1 = 2 * (g1 * h2 + g2 * h1)

    R = (
        (g1**2 + h1**2) * exp(2 * alpha1)
        + (g2**2 + h2**2) * exp(-2 * alpha1)
        + A * cos(2 * gamma1)
        + B * sin(2 * gamma1)
    )
    R /= (
        exp(2 * alpha1)
        + (g1**2 + h1**2) * (g2**2 + h2**2) * exp(-2 * alpha1)
        + C1 * cos(2 * gamma1)
        + D1 * sin(2 * gamma1)
    )
    R = R.real

    T = (n2 / n0) * ((1 + g1) ** 2 + h1**2) * ((1 + g2) ** 2 + h2**2)
    T /= (
        exp(2 * alpha1)
        + (g1**2 + h1**2) * (g2**2 + h2**2) * exp(-2 * alpha1)
        + C1 * cos(2 * gamma1)
        + D1 * sin(2 * gamma1)
    )
    T = T.real

    return R, T, 1 - R - T
