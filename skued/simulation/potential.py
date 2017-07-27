# -*- coding: utf-8 -*-
"""
Electrostatic potential simulation
==================================
"""
from functools import partial
from .. import chunked, minimum_image_distance, repeated_array
from scipy.special import k0 as bessel
import numpy as np
from numpy import pi

m = 9.109*10**(-31)     #in kg
a0 = 0.5291             #in Angs
e = 14.4                #Volt*Angstrom

def electrostatic(crystal, x, y, z):
    """
    Electrostatic potential from a crystal calculated on a real-space mesh, 
    assuming an infinite crystal.

    Parameters
    ----------
    crystal : skued.Crystal
        
    x, y, z : `~numpy.ndarray`
        Real space coordinates mesh. 
    
    Returns
    -------
    potential : `~numpy.ndarray`, dtype float
        Linear superposition of atomic potential [V*Angs]
    
    See also
    --------
    pelectrostatic: projected electrostatic potential of an infinite crystal.
    """
    # TODO: split xx and yy into smalled non-repeating unit
    # TODO: multicore
    
    potential = np.zeros_like(x, dtype = np.float)
    r = np.zeros_like(x, dtype = np.float)
    for atom in crystal:
        ax, ay, az = atom.xyz(crystal)
        r[:] = minimum_image_distance(x - ax, y - ay, z - az, 
                                        lattice = crystal.lattice_vectors)
        potential += atom.potential(r)
    
    # Due to sampling, x,y, and z might pass through the center of atoms
    # Replace np.inf by the next largest value
    m = potential[np.isfinite(potential)].max()
    potential[np.isinf(potential)] = m
    return potential

def pelectrostatic(crystal, x, y, bounds = None):
    """
    Projected electrostatic potential from a crystal calculated on a real-space mesh, 
    assuming an infinite crystal in x and y. Projection axis is defined as the z-axis. 
    To project the potential along a different axis, the crystal can be rotated with ``Crystal.transform``.

    Parameters
    ----------
    crystal : skued.Crystal
        
    x, y:  `~numpy.ndarray`
        Real-space coordinates. 
    bounds : iterable or None, optional
        Bounds of atom inclusion. Atoms with real-space z-position outside [ min(bounds), max(bounds) )
        are not counted in the computation.
    
    Returns
    -------
    potential : `~numpy.ndarray`, dtype float
        Linear superposition of electrostatic potential [V*Angs]
    
    See also
    --------
    electrostatic: three-dimensional electrostatic potential of an infinite crystal.
    """
    # TODO: split xx and yy into smalled non-repeating unit
    #       np.unique(np.mod(xx, per_x))
    # TODO: multicore

    if bounds:
        min_z, max_z = min(bounds), max(bounds)
        atoms = (atom for atom in iter(crystal) if min_z <= atom.xyz(crystal)[2] < max_z)
    else:
        atoms = iter(crystal)

    potential = np.zeros_like(x, dtype = np.float)
    z = np.zeros_like(x)
    for atom in atoms:
        xa, ya, _ = atom.xyz(crystal)
        r = minimum_image_distance(x - xa, y - ya, z, lattice = np.array(crystal.lattice_vectors)).reshape(-1,1)
        potential += np.sum( 2*atom._a*bessel(2*pi*r*np.sqrt(atom._b)) + (atom._c/atom._d) * np.exp( -(r*pi)**2 / atom._d), axis = -1).reshape(x.shape)
    potential *= 2 * a0 * e * (pi**2)
    
    # Due to sampling, x,y, and z might pass through the center of atoms
    # Replace n.inf by the next largest value
    potential[np.isinf(potential)] = np.nan
    potential[np.isnan(potential)] = np.nanmax(potential)
    return potential