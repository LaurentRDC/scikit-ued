# -*- coding: utf-8 -*-
"""
Structure Factor calculation
"""
from collections import Iterable
from functools import lru_cache
from itertools import count, product, takewhile
from math import sqrt

import numpy as np
from numpy.linalg import norm
from npstreams import primed

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
    Gx, Gy, Gz = crystal.scattering_vector(h, k, l)
    nG = np.sqrt(Gx ** 2 + Gy ** 2 + Gz ** 2)

    # Separating the structure factor into sine and cosine parts avoids adding
    # complex arrays together. About 3x speedup vs. using complex exponentials
    SFsin, SFcos = (
        np.zeros(shape=nG.shape, dtype=np.float),
        np.zeros(shape=nG.shape, dtype=np.float),
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

# Generalized hypotenuse
def _hypot(*args):
    return sqrt(sum(map(lambda i: i**2, args)))

# We prime this generator so that pre-checks are evaluated
# For example, no negative bounds.
# This is better for quick feedback while writing a script.
@primed
def bounded_reflections(crystal, nG):
    """
    Generates reflections (hkl) with norm(G) <= nG
    
    Parameters
    ----------
    crystal : Crystal
        Crystal instance
    nG : float
        Maximal scattering vector norm. By our convention, :math:`G = 4 \pi s`.
    
    Yields
    ------
    reflection : 3-tuple of ints
        Miller indices of a bounded reflection.
    """
    if nG < 0:
        raise ValueError("Bound {} is negative.".format(nG))
    
    yield
    # Determine the maximum index such that (i00) family is still within data limits
    # This provides a (large) upper bound so that we are sure that the overall filtering will terminate
    bounded = lambda i: any(
        [
            norm(crystal.scattering_vector(i, 0, 0)) <= nG,
            norm(crystal.scattering_vector(0, i, 0)) <= nG,
            norm(crystal.scattering_vector(0, 0, i)) <= nG,
        ]
    )
    max_index = max(takewhile(bounded, count(0)))
    extent = range(-max_index, max_index + 1)
    refls = product(extent, repeat=3)

    # The above bound was only a first pass. We can refine further
    in_bounds = lambda refl: _hypot(*crystal.scattering_vector(*refl)) <= nG
    yield from filter(in_bounds, refls)