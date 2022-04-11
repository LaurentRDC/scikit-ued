# -*- coding: utf-8 -*-
""" Electron properties """

import numpy as np
from scipy.constants import Planck, electron_mass, elementary_charge, speed_of_light


def lorentz(keV):
    """
    Relativistic factor :math:`\\gamma`, defined as :math:`\\gamma = \\frac{1}{\\sqrt{1 - v^2/c^2}}`

    Parameters
    ----------
    keV : array_like or float
        Electron energy [keV].

    Returns
    -------
    out : array_like or float

    References
    ----------
    .. Kirkland 2010 Eq. 2.2
    """
    return 1 + (elementary_charge * keV * 1e3) / (electron_mass * speed_of_light**2)


def electron_wavelength(keV):
    """
    Relativistic wavelength :math:`\\lambda` of an accelerated electron.

    .. math::

        \\lambda = \\frac{h}{\\sqrt{2 m_e e V}}\\gamma

    where :math:`\\gamma` is the relativistic Lorentz factor.

    Parameters
    ----------
    keV : array_like or float
        Electron energy [keV].

    Returns
    -------
    out : array_like or float
        Electron wavelength [:math:`Å^{-1}`]

    References
    ----------
    .. Kirkland 2010 Eq. 2.5
    """
    eV = elementary_charge * keV * 1e3
    wavelength_meters = (
        Planck
        * speed_of_light
        / np.sqrt(eV * (2 * electron_mass * speed_of_light**2 + eV))
    )
    return wavelength_meters * 1e10  # wavelength in angstroms


def electron_velocity(keV):
    """
    Relativistic velocity :math:`v_e` of an accelerated electron.

    .. math::

        \\frac{v_e}{c} = \\sqrt{1 - \\frac{m_0 c^2}{m_0 c^2 + e V}}

    Parameters
    ----------
    keV : array_like or float
        Electron energy [keV].

    Returns
    -------
    out : array_like or float
        Electron velocity [:math:`Å / s`]

    References
    ----------
    .. Kirkland 2010 Eq. 2.3
    """
    eV = elementary_charge * keV * 1e3
    m0c2 = electron_mass * speed_of_light**2
    v_over_c = np.sqrt(eV * (eV + 2 * m0c2)) / (m0c2 + eV)
    return (speed_of_light * v_over_c) * 1e10  # speed in Angstroms


def interaction_parameter(keV):
    """
    Interaction parameter from relativistic electron wavelength.

    Parameters
    ----------
    keV : array_like or float
        Electron energy [keV].

    Returns
    -------
    out : array_like or float
        Interaction parameter [:math:`rad/(V Å)`]

    References
    ----------
    .. Kirkland 2010 Eq. 5.6
    """
    l = electron_wavelength(keV)
    V = keV * 1e3

    return (
        (2 * np.pi)
        / (electron_wavelength(keV) * V)
        * (electron_mass * speed_of_light**2 + elementary_charge * V)
        / (2 * electron_mass * speed_of_light**2 + elementary_charge * V)
    )
