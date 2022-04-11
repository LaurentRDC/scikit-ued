# -*- coding: utf-8 -*-
"""
Structure Factor calculation
"""

import numpy as np

from crystals.affine import change_basis_mesh

from .form_factors import affe


def structure_factor(crystal, h, k, l, normalized=False):
    """
    Computation of the static structure factor for electron diffraction.

    Parameters
    ----------
    crystal : Crystal
        Crystal instance
    h, k, l : array_likes or floats
        Miller indices. Can be given in a few different formats:

        * floats : returns structure factor computed for a single scattering vector

        * 3 coordinate ndarrays, shapes (L,M,N) : returns structure factor computed over all coordinate space

    normalized : bool, optional
        If True, the normalized structure factor :math`E` is returned. This is the statis structure
        factor normalized by the sum of form factors squared.

    Returns
    -------
    sf : ndarray, dtype complex
        Output is the same shape as input G[0]. Takes into account
        the Debye-Waller effect.
    """
    # Distribute input
    # This works whether G is a list of 3 numbers, a ndarray shape(3,) or
    # a list of meshgrid arrays.
    h, k, l = np.atleast_1d(h, k, l)
    Gx, Gy, Gz = change_basis_mesh(
        h, k, l, basis1=crystal.reciprocal_vectors, basis2=np.eye(3)
    )
    nG = np.sqrt(Gx**2 + Gy**2 + Gz**2)

    # Separating the structure factor into sine and cosine parts avoids adding
    # complex arrays together. About 3x speedup vs. using complex exponentials
    SFsin, SFcos = (
        np.zeros(shape=nG.shape, dtype=float),
        np.zeros(shape=nG.shape, dtype=float),
    )

    # Pre-allocation of form factors gives huge speedups
    atomff_dict = dict()
    for atom in crystal:  # TODO: implement in parallel?

        if atom.element not in atomff_dict:
            atomff_dict[atom.element] = affe(atom, nG)

        x, y, z = atom.coords_cartesian
        arg = x * Gx + y * Gy + z * Gz
        # TODO: debye waller factor based on displacement
        atomff = atomff_dict[atom.element]
        SFsin += atomff * np.sin(arg)
        SFcos += atomff * np.cos(arg)

    SF = SFcos + 1j * SFsin

    if normalized:
        SF /= np.sqrt(sum(atomff_dict[atom.element] ** 2 for atom in crystal))

    return SF
