# -*- coding: utf-8 -*-
"""
Electron scattering form factors
"""

from pathlib import Path

import numpy as np
from yaml import load
from math import sin, cos

from crystals import Element, Atom

from .scattering_params import scattering_params

# For aspherical e form factors
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

DATADIR = Path(__file__).parent / "data"

with open(DATADIR / "aspherical.yaml") as f:
    aspherical_ff = load(f, Loader=Loader)


def affe(atom, nG):
    """
    Atomic form factors for electrons, for neutral atoms.

    Parameters
    ----------
    atom : skued.Atom, int, or str
        Atomic number, atomic symbol, or Atom instance.
        Atomic symbols are expected to be properly capitalized, e.g. ``As`` or ``W``.
    nG : array_like
        Scattering vector norm, in units of Angstroms:math:`^{-1}`. (:math:`|G| = 4 \\pi s`).

    Returns
    -------
    eff : `~numpy.ndarray`, dtype float
        Atomic form factor for electron scattering.

    Raises
    ------
    ValueError : scattering information is not available, for example if the atomic number is larger than 103
    """
    if isinstance(atom, (int, str)):
        atom = Element(atom)
    atomic_number = atom.atomic_number

    try:
        _, a1, b1, a2, b2, a3, b3, c1, d1, c2, d2, c3, d3 = scattering_params[
            atomic_number
        ]
    except KeyError:
        raise ValueError(
            f"Scattering information for element Z={atomic_number} is unavailable."
        )

    # Parametrization of form factors is done in terms of q = 2 s = 2 pi |G|
    q = nG / (2 * np.pi)
    q2 = np.square(q)
    sum1 = a1 / (q2 + b1) + a2 / (q2 + b2) + a3 / (q2 + b3)
    sum2 = c1 * np.exp(-d1 * q2) + c2 * np.exp(-d2 * q2) + c3 * np.exp(-d3 * q2)
    return sum1 + sum2


def aspherical_affe(atom, s):
    """
    Aspherical atomic form factors for electron scattering.
    Only atoms lighter than Xe (and including) are supported (Z <= 54).

    Parameters
    ----------
    atom : crystals.Atom, int, or str
        Atomic number, atomic symbol, or Atom instance.
        Atomic symbols are expected to be properly capitalized, e.g. ``As`` or ``W``.
        Note that if `atom` is provided as an int or str, the ground state electronic structure
        is used. If an `crystals.Atom` object is provided, then its electronic structure will
        be used to construct the electron form factor.
    s : array_like
        Scattering vector norm

    Returns
    -------
    eff : `~numpy.ndarray`, dtype float
        Atomic form factor for electron scattering.

    Raises
    ------
    ValueError : scattering information is not available, for example if the atomic number is larger than 54

    References
    ----------
    .. [#] Jin-Cheng Zheng, Lijun Wu and Yimei Zhu. "Aspherical electron scattering factors and their
           parameterizations for elements from H to Xe" (2009). J. Appl. Cryst. vol 42, pp. 1043 - 1053.
    """
    if isinstance(atom, (int, str)):
        atom = Atom(atom)
    element = atom.element

    return _affe_parametrization(s, aspherical_ff[element]["total"])


def _affe_parametrization(s, d):
    """Reconstruct affe parametrization.

    Parameters
    ----------
    s : array_like
        Scattering vector norm
    d : dict[str, array_like]
        Dictionary containing the "a" and "b" parameters.
    """
    s2 = np.square(np.asfarray(s))
    result = np.zeros_like(s2)
    for a, b in zip(d["a"], d["b"]):
        result += a * np.exp(-b * s2)
    return result


def _affe_p(element, s):
    """
    Atomic form factor for p orbitals. Effectively, this is from the reference, Eq. 11

    .. [#] Jin-Cheng Zheng, Lijun Wu and Yimei Zhu. "Aspherical electron scattering factors and their
           parameterizations for elements from H to Xe" (2009). J. Appl. Cryst. vol 42, pp. 1043 - 1053.
    """
    affe_p_sph = _affe_parametrization(
        s, aspherical_ff[element]["p0"]
    )  # Numerical parametrization of Eq 12
    affe_p_p1 = _affe_parametrization(
        s, aspherical_ff[element]["p1"]
    )  # Numerical parametrization of Eq 13

    # Angle between the electron beam and the z-axis of the orbital
    # TODO: How do we find this?
    beta = 0  # radians
    return (3 / 2) * (sin(beta) ** 2) * affe_p_sph + (
        cos(beta) ** 2 - (1 / 2) * sin(beta) ** 2
    ) * affe_p_p1
