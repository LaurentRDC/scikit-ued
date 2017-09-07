# -*- coding: utf-8 -*-
"""
Polycrystalline diffraction pattern simulation
==============================================
"""
from ..voigt import pseudo_voigt
from .structure_factors import structure_factor, bounded_reflections
import numpy as np

def powdersim(crystal, scattering_length, fwhm_g = 0.01, fwhm_l = 0.02, **kwargs):
    """
    Simulates polycrystalline diffraction pattern.

    Parameters
    ----------
    crystal : `skued.structure.Crystal`
        Crystal from which to diffract.
    scattering_length : `~numpy.ndarray`, shape (N,)
        Range of scattering length over which to compute the diffraction pattern [2pi/Angs].
    fwhm_g, fwhm_l : float, optional
            Full-width at half-max of the Gaussian and Lorentzian parts of the Voigt profile.
            See `skued.pseudo_voigt` for more details.

    Returns
    -------
    pattern : `~numpy.ndarray`, shape (N,)
        Diffraction pattern
    """
    h, k, l = bounded_reflections(crystal, nG = 4*np.pi*scattering_length.max())
    Gx, Gy, Gz = crystal.scattering_vector(h, k, l)
    scatt_length = np.sqrt(Gx**2 + Gy**2 + Gz**2)/(4*np.pi)
    intensities = np.absolute(structure_factor(crystal, h, k, l))**2

    psf = pseudo_voigt(scattering_length, center = np.mean(scattering_length), fwhm_g = fwhm_g, fwhm_l = fwhm_l)

    pattern = np.zeros_like(scattering_length)
    for s, i in zip(scatt_length, intensities):
        pattern += i * pseudo_voigt(scattering_length, s, fwhm_g, fwhm_l)

    return pattern