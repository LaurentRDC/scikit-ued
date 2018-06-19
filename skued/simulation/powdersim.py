# -*- coding: utf-8 -*-
"""
Polycrystalline diffraction pattern simulation
==============================================
"""
import numpy as np

from ..voigt import pseudo_voigt
from .structure_factors import bounded_reflections, structure_factor


def powdersim(crystal, q, fwhm_g = 0.03, fwhm_l = 0.06, **kwargs):
    """
    Simulates polycrystalline diffraction pattern.

    Parameters
    ----------
    crystal : `skued.structure.Crystal`
        Crystal from which to diffract.
    q : `~numpy.ndarray`, shape (N,)
        Range of scattering vector norm over which to compute the diffraction pattern [1/Angs].
    fwhm_g, fwhm_l : float, optional
        Full-width at half-max of the Gaussian and Lorentzian parts of the Voigt profile.
        See `skued.pseudo_voigt` for more details.

    Returns
    -------
    pattern : `~numpy.ndarray`, shape (N,)
        Diffraction pattern
    """
    h, k, l = bounded_reflections(crystal, nG = q.max())
    Gx, Gy, Gz = crystal.scattering_vector(h, k, l)
    qs = np.sqrt(Gx**2 + Gy**2 + Gz**2)
    intensities = np.absolute(structure_factor(crystal, h, k, l))**2

    pattern = np.zeros_like(q)
    for qi, i in zip(qs, intensities):
        pattern += i * pseudo_voigt(q, qi, fwhm_g, fwhm_l)

    return pattern
