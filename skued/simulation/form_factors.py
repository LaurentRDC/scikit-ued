# -*- coding: utf-8 -*-
"""
Electron scattering form factors
"""

from pathlib import Path

import numpy as np
from yaml import load

from ..structure import NUM_TO_ELEM, ELEM_TO_NUM
from .scattering_params import scattering_params

# For aspherical e form factors
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

DATADIR = Path(__file__).parent / "data"

with open(DATADIR / "aspherical.yaml") as f:
    aspherical_ff = load(f, Loader=Loader)


def aspherical_affe(atom, s):
    """ 
    Aspherical atomic form factors for electron scattering. 
    Only atoms lighter than Xe (and including) are supported (Z <= 54).

    Parameters
    ----------
    atom : skued.Atom, int, or str
        Atomic number, atomic symbol, or Atom instance.
        Atomic symbols are expected to be properly capitalized, e.g. ``As`` or ``W``.
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
    if isinstance(atom, int):
        element = NUM_TO_ELEM[atom]
    elif isinstance(atom, str):
        element = atom
    else:
        element = atom.element

    params_a = aspherical_ff[element]["total"]["a"]
    params_b = aspherical_ff[element]["total"]["b"]

    s2 = np.square(np.asfarray(s))

    result = np.zeros_like(s2)
    for a, b in zip(params_a, params_b):
        result += a * np.exp(-b * s2)
    return result


def affe(atom, nG):
    """
    Atomic form factors for electrons, for neutral atoms. 

    Parameters
    ----------
    atom : skued.Atom, int, or str
        Atomic number, atomic symbol, or Atom instance.
        Atomic symbols are expected to be properly capitalized, e.g. ``As`` or ``W``.
    nG : array_like
        Scattering vector norm, in units of Angstroms:math:`^{-1}`. (:math:`|G| = 4 \pi s`). 
    
    Returns
    -------
    eff : `~numpy.ndarray`, dtype float
        Atomic form factor for electron scattering.

    Raises
    ------
    ValueError : scattering information is not available, for example if the atomic number is larger than 103
    """
    if isinstance(atom, int):
        atomic_number = atom
    elif isinstance(atom, str):
        atomic_number = ELEM_TO_NUM[atom]
    else:
        atomic_number = atom.atomic_number

    try:
        _, a1, b1, a2, b2, a3, b3, c1, d1, c2, d2, c3, d3 = scattering_params[
            atomic_number
        ]
    except KeyError:
        raise ValueError(
            "Scattering information for element Z={} is unavailable.".format(
                atomic_number
            )
        )

    # Parametrization of form factors is done in terms of q = 2 s = 2 pi |G|
    q = nG / (2 * np.pi)
    q2 = np.square(q)
    sum1 = a1 / (q2 + b1) + a2 / (q2 + b2) + a3 / (q2 + b3)
    sum2 = c1 * np.exp(-d1 * q2) + c2 * np.exp(-d2 * q2) + c3 * np.exp(-d3 * q2)
    return sum1 + sum2
