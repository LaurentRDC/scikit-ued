# -*- coding: utf-8 -*-
"""
Structure Factor calculation
"""
from collections import Iterable
from itertools import takewhile, count, product

import numpy as np
from numpy.linalg import norm

from .. import change_basis_mesh
from ..structure import ELEM_TO_NUM
from .scattering_params import scattering_params

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
        _, a1, b1, a2, b2, a3, b3, c1, d1, c2, d2, c3, d3 = scattering_params[atomic_number]
    except KeyError:
        raise ValueError('Scattering information for element Z={} is unavailable.'.format(atomic_number))
    
    # Parametrization of form factors is done in terms of q = 2 s = 2 pi |G|
    q = nG / (2*np.pi)
    q2 = np.square(q)
    sum1 = a1/(q2 + b1) + a2/(q2 + b2) + a3/(q2 + b3)
    sum2 = c1 * np.exp(-d1 * q2) + c2 * np.exp(-d2 * q2) + c3 * np.exp(-d3 * q2)
    return sum1 + sum2

def structure_factor(crystal, h, k, l, normalized = False):
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
    nG = np.sqrt(Gx**2 + Gy**2 + Gz**2)
    
    # Separating the structure factor into sine and cosine parts avoids adding
    # complex arrays together. About 3x speedup vs. using complex exponentials
    SFsin, SFcos = np.zeros(shape = nG.shape, dtype = np.float), np.zeros(shape = nG.shape, dtype = np.float)

    # Pre-allocation of form factors gives huge speedups
    dwf = np.empty_like(SFsin) 	# debye-waller factor
    atomff_dict = dict()
    for atom in crystal: #TODO: implement in parallel?

        if atom.element not in atomff_dict:
            atomff_dict[atom.element] = affe(atom, nG)

        x, y, z = atom.xyz(crystal)
        arg = x*Gx + y*Gy + z*Gz
        atom.debye_waller_factor((Gx, Gy, Gz), out = dwf)
        atomff = atomff_dict[atom.element]
        SFsin += atomff * dwf * np.sin(arg)
        SFcos += atomff * dwf * np.cos(arg)
    
    SF = SFcos + 1j*SFsin

    if normalized:
        SF /= np.sqrt(sum(atomff_dict[atom.element]**2 for atom in crystal))
    
    return SF

def bounded_reflections(crystal, nG):
    """
    Returns iterable of reflections (hkl) with norm(G) < nG
    
    Parameters
    ----------
    crystal : Crystal
        Crystal instance
    nG : float
        Maximal scattering vector norm. By our convention, :math:`G = 4 \pi s`.
    
    Returns
    -------
    h, k, l : ndarrays, shapes (N,), dtype int
    """
    if nG < 0:
        raise ValueError('Bound {} is negative.'.format(nG))
    
    # Determine the maximum index such that (i00) family is still within data limits
    #TODO: cache results based on max_index?
    bounded = lambda i : any([norm(crystal.scattering_vector(i,0,0)) <= nG, 
                                norm(crystal.scattering_vector(0,i,0)) <= nG, 
                                norm(crystal.scattering_vector(0,0,i)) <= nG])
    max_index = max(takewhile(bounded, count(0)))
    extent = range(-max_index, max_index + 1)
    h, k, l = np.split(np.array(list(product(extent, extent, extent)), dtype = np.int), 3, axis = -1)
    h, k, l = h.ravel(), k.ravel(), l.ravel()

    # we only have an upper bound on possible reflections
    # Let's filter down
    Gx, Gy, Gz = crystal.scattering_vector(h, k, l)
    norm_G = np.sqrt(Gx**2 + Gy**2 + Gz**2)
    in_bound = norm_G <= nG
    return h.compress(in_bound), k.compress(in_bound), l.compress(in_bound)