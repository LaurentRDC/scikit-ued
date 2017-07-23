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
    pelectrostatic
        Projected electrostatic potential from a crystal
    """
    # TODO: split xx and yy into smalled non-repeating unit
    # TODO: multicore
    # TODO: pre-load scattering params _a, _b, _c, and _d into an appropriate shape
    
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

def pelectrostatic(crystal, xx, yy, bounds = None):
    """
    Projected electrostatic potential from a crystal calculated on a real-space mesh, 
    assuming an infinite crystal. Projection axis is defined as the z-axis. To project the potential
    along a different axis, the crystal can be rotated with ``Crystal.transform``.

    Parameters
    ----------
    crystal : skued.Crystal
        
    xx, yy:  `~numpy.ndarray`
        Real space coordinates mesh. 
    
    Returns
    -------
    potential : `~numpy.ndarray`, dtype float
        Linear superposition of atomic potential [V*Angs]
    """
    # TODO: split xx and yy into smalled non-repeating unit
    # TODO: multicore
    # TODO: pre-load scattering params _a, _b, _c, and _d into an appropriate shape

    if bounds:
        min_z, max_z = min(bounds), max(bounds)
        atoms = (atom for atom in iter(crystal) if min_z <= atom.xyz(crystal)[2] < max_z)
    else:
        atoms = iter(crystal)

    potential = np.zeros_like(xx)
    zz = np.zeros_like(xx)
    for atom in atoms:
        xa, ya, _ = atom.xyz(crystal)
        r = minimum_image_distance(xx - xa, yy - ya, zz, lattice = np.array(crystal.lattice_vectors)).reshape(-1,1)
        potential += np.sum( 2*atom._a*bessel(2*pi*r*np.sqrt(atom._b)) + (atom._c/atom._d) * np.exp( -(r*pi)**2 / atom._d), axis = -1).reshape(xx.shape)
    potential *= 2 * a0 * e * (pi**2)
    
    # Due to sampling, x,y, and z might pass through the center of atoms
    # Replace n.inf by the next largest value
    potential[np.isinf(potential)] = np.nan
    potential[np.isnan(potential)] = np.nanmax(potential)
    return potential

def _proj_elec_atm(atom, xx, yy, lattice):
    potential = np.zeros_like(xx, dtype = np.float)[:,:,None]
    xa, ya, _ = tuple(atom.xyz(lattice))
    r = np.zeros_like(xx)[:,:,None, None]
    r[:,:, 0, 0] = minimum_image_distance(xx - xa, yy - ya, np.zeros_like(xx), lattice)
    potential = np.sum( 2*atom._a*bessel(2*pi*r*np.sqrt(atom._b)) + (atom._c/atom._d) * np.exp( -(r*pi)**2 / atom._d), axis = -1)
    return 2*a0*e*(pi**2)*np.squeeze(potential)   