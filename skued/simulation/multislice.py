# -*- coding: utf-8 -*-
"""
Diffraction simulation
======================
"""

from functools import partial
from itertools import chain, repeat

import numpy as np

from .potential import pelectrostatic
from .fourier import fft2freq, limit_bandwidth, next_fast_shape
from .. import electron_wavelength, interaction_parameter

def multislice(crystal, thickness = None, energy = 90, shape = (512, 512), resolution = 1):
    """
    Propagate a plane wave through a crystal using the multislice algorithm.

    Parameters
    ----------
    crystal : skued.Crystal
        Crystal from which to scatter.
    thickness : float or None, optional
        Thickness of the sample in angstroms. The thickness will be rounded to the nearest
        unit-cell thickness. If None (default), the sample will be propagated through a single
        unit cell, equivalent to the weak-phase approximation.
    energy : float, optional
        Electron energy [keV]
    shape : 2-tuple of ints, optional
        Shape of the output, which will be rounded up to the nearest fast FFT shape.
        The simulation supercell dimensions will be determined based on this output shape.
    resolution : float, optional
        Desired simulation resolution [:math:`Å^{-1}`].

    Returns
    -------
    exit_wave : `~numpy.ndarray`
        Scattered wave.
    
    Raises
    ------
    RuntimeError : if propagation intensity decreases by more than 10%.
    """
    shape = next_fast_shape(shape)
    period_x, period_y, period_z = crystal.periodicity

    # Determine the real-space resolution
    # Note that we require that the final extent be a multiple of the periodicity
    dx, dy = resolution/2, resolution/2
    x_extent, y_extent = dx * shape[0], dy * shape[1]   # Might not be periodic
    x_extent = round(x_extent / period_x) * period_x
    y_extent = round(y_extent / period_y) * period_y
    extent = 10*period_x*period_y
    
    # Determine slice thickness as a divider of the unitcell period in z
    # We try to have the slice thickness be close to 1 angstroms
    if thickness is None:
        n_slices = 1    # weak-phase approximation
    else:
        n_slices = int(period_z)

    if thickness is None:
        thickness = period_z

    xx, yy = np.meshgrid(np.linspace(0, extent, num = shape[0]), 
                         np.linspace(0, extent, num = shape[1]), indexing = 'xy')
    kx, ky = fft2freq(xx, yy, indexing = 'xy')
    k = np.hypot(kx, ky)

    # Calculate transmission functions through the slices
    # TODO: is it enough to only bandwidth-limit the potential slices or must this be
    #       done for transmission functions as well?
    slices = [limit_bandwidth(np.exp(1j * interaction_parameter(energy) * p), k, (1/2)*k.max()) 
              for p in potential_slices(crystal, xx, yy, n_slices = n_slices)]
    
    wavefunction = np.ones_like(xx, dtype = np.complex)
    initial_intensity = np.sum(np.abs(wavefunction)**2)
    propagator = np.exp(1j * np.pi * electron_wavelength(energy) * (period_z / n_slices) * k**2)

    # Cycle through the slices until full thickness
    n_cells = int(thickness / period_z)     # number of unit cells in transverse direction to go through
    cycled_slices = chain.from_iterable(repeat(slices, n_cells))
    for transmission in cycled_slices:
        wavefunction = np.fft.ifft2( propagator * np.fft.fft2(wavefunction * transmission))

        relative_intensity = np.sum(np.abs(wavefunction)**2)/initial_intensity
        if relative_intensity < 0.9:
            raise RuntimeError('Simulation did not converge; the wavefunction intensity has decreased by more than 10%.')

    # Apply objective lens transfer function
    # TODO: allow changing parameters
    obj_tf = objective_TF(k, energy)
    wavefunction = np.fft.ifft2(np.fft.fft2(wavefunction) * obj_tf)

    return wavefunction

def objective_TF(k, energy, defocus = 0, spherical_aberration = 0):
    """
    Objective lens transfer function.

    Parameters
    ----------
    k : `~numppy.ndarray`, ndim 2
        Wavevector norm
    energy : float, optional
        Electron energy [keV]
    defocus : float
        Lens defocus parameter [:math:`Å^{-1}`].
    spherical_aberration : float
        Spherical aberration parameter.
    
    Returns
    -------
    transfer_function : `~numppy.ndarray`, ndim 2
        Transfer function of the objective lens.
    """
    # TODO: aperture function, i.e. max semiangle
    aperture_TF = np.ones_like(k)

    wavelength = electron_wavelength(energy)
    k2 = np.square(k)
    chi = np.pi * wavelength * k2 * (0.5 * spherical_aberration * wavelength**2 * k2 - defocus)
    return np.exp(1j*chi) * aperture_TF

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