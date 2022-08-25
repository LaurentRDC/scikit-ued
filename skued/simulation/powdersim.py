# -*- coding: utf-8 -*-
"""
Polycrystalline diffraction pattern simulation
==============================================
"""
import numpy as np

from crystals.affine import change_basis_mesh

from ..voigt import pseudo_voigt
from .structure_factors import structure_factor


def powdersim(crystal, q, fwhm_g=0.03, fwhm_l=0.06, **kwargs):
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
    refls = np.vstack(tuple(crystal.bounded_reflections(q.max())))
    h, k, l = np.hsplit(refls, 3)
    Gx, Gy, Gz = change_basis_mesh(
        h, k, l, basis1=crystal.reciprocal_vectors, basis2=np.eye(3)
    )
    qs = np.sqrt(Gx**2 + Gy**2 + Gz**2)
    intensities = np.absolute(structure_factor(crystal, h, k, l)) ** 2

    pattern = np.zeros_like(q)
    for qi, i in zip(qs, intensities):
        pattern += i * pseudo_voigt(q, qi, fwhm_g, fwhm_l)

    return pattern
