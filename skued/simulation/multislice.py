# -*- coding: utf-8 -*-
"""
Diffraction simulation
======================
"""

from functools import partial
from itertools import chain, repeat

import numpy as np

from .potential import pelectrostatic
from .fourier import fft2freq, limit_bandwidth
from .. import electron_wavelength, interaction_parameter

def multislice(crystal, thickness = None, energy = 90):
    """
    Propagate a plane wave through a crystal using the multislice algorithm.

    Parameters
    ----------
    crystal : skued.Crystal
        Crystal from which to scatter.
    thickness : float or None, optional
        Thickness of the sample in angstroms. If None (default), the sample will be propagated through a single
        unit cell, equivalent to the weak-phase approximation.
    energy : float, optional
        Electron energy [keV]
    
    Returns
    -------
    exit_wave : `~numpy.ndarray`
        Scattered wave.
    """
    # TODO: allow adjustable unit cell slicing thickness (i.e. n_slices)
    # TODO: allow output shape
    if thickness is None:
        n_slices = 1    # weak-phase approximation
    else:
        n_slices = 10
    shape = (512, 512)

    # Determine slice thickness as a divider of the unitcell period in z
    # We try to have the slice thickness be close to 2 angstroms
    # TODO: dynamically determine the perfect thickness?
    period_x, period_y, period_z = crystal.periodicity
    extent = 10*period_x*period_y

    if thickness is None:
        thickness = period_z

    xx, yy = np.meshgrid(np.linspace(0, extent, num = shape[0]), 
                         np.linspace(0, extent, num = shape[1]), indexing = 'xy')
    kx, ky = fft2freq(xx, yy, indexing = 'xy')
    k = np.hypot(kx, ky)

    # Calculate transmission functions through the slices
    # TODO: is it enough to only bandwidth-limit the potential slices or must this be
    #       done for transmission functions as well?
    slices = [limit_bandwidth(p, k, (2/3)*k.max()) for p in potential_slices(crystal, xx, yy, n_slices = n_slices)]
    
    wavefunction = np.ones_like(xx, dtype = np.complex)
    propagator = np.exp(1j * np.pi * electron_wavelength(energy) * (period_z / n_slices) * k**2)

    n_steps = int(thickness / period_z)             # number of iterations
    for potential_slice in ncycles(slices, n_steps):
        transmission = np.exp(1j * interaction_parameter(energy) * potential_slice)
        transmission[:] = limit_bandwidth(transmission, k, (2/3)*k.max())
        wavefunction = np.fft.ifft2( propagator * np.fft.fft2(wavefunction * transmission))

    return wavefunction

def ncycles(iterable, n):
    """Returns the sequence elements n times"""
    return chain.from_iterable(repeat(tuple(iterable), n))

def potential_slices(crystal, X, Y, n_slices):
    """
    Generate electrostatic potential slices for a unit cell.

    Parameters
    ----------
    crystal : skued.Crystal
        
    x, y:  `~numpy.ndarray`
        Real-space coordinates.
    n_slices : int
        Number of slices to computed inside the unit cell.
    
    Yields
    ------
    potential : `~numpy.ndarray`
    """
    # Compute the reciprocal space coordinates for better bandwidth limitation
    _, _, period_z = crystal.periodicity
    slice_thickness = period_z / n_slices

    span = (0, slice_thickness)
    for slice_index in range(n_slices):
        slice_start = slice_index * slice_thickness
        potential = pelectrostatic(crystal, X, Y, bounds = (slice_start, slice_start + slice_thickness))

        yield potential