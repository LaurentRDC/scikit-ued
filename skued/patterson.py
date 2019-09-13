# -*- coding: utf-8 -*-
"""
Patterson pair-correlation function for polycrystalline diffraction patterns
============================================================================
"""

import numpy as np

from .simulation import affe


def patterson(q, I, crystal, radii):
    """
    Computation of the patterson pair-distribution function from azimuthal diffraction pattern.

    Parameters
    ----------
    q : ndarray, shape (N,)
        Scattering vector associated with the diffracted intensity `I` [1/$\AA$].
    I : ndarray, shape (N,)
        Diffracted intensity.
    crystal : crystals.Crystal
        Crystal structure associated with this data.
    radii : ndarray, shape (M,)
        Array over which the pair-correlation function is calculated [$\AA$].
    
    Returns
    -------
    patterson : ndarray, shape (M,)
        Patterson pair-distribution function calculated over the array `radii`.
    
    References
    ----------
    [1] Fultz and Howe, Transmission Electron Microscopy and DIffractometry of Materials, equation 9.144
    """
    q, I, radii = tuple(map(np.asfarray, [q, I, radii]))

    normalization = sum(np.square(affe(atm, q)) for atm in crystal)
    reduced_intensity = I / normalization

    # We employ the outer product to avoid loops
    # integral over q is done over axis 1
    # Therefore, we extend the reduced intensity over axis 1
    rr, GG = np.meshgrid(radii, q)
    extended_reduced_intensity = np.outer(reduced_intensity, np.ones_like(radii))

    dq = np.mean(np.diff(q))
    return (
        np.sum((GG / rr) * np.sin(rr * GG) * (extended_reduced_intensity - 1), axis=0)
        * dq
    )