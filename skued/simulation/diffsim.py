# -*- coding: utf-8 -*-
"""
Diffraction simulation
======================
"""
from itertools import chain, repeat
import numpy as np
from math import ceil, atan
from scipy.fftpack import fftfreq, fftshift
from scipy.interpolate import RegularGridInterpolator
from warnings import warn

from . import pelectrostatic
from .. import interaction_parameter, electron_wavelength

FFTOPS = {}
try:
    from pyfftw.interfaces.scipy_fftpack import fft2, ifft2
    from pyfftw.interfaces.cache import enable, set_keepalive_time
    FFTOPS['threads'] = 2
    enable()
    set_keepalive_time(1)
except ImportError:
    from scipy.fftpack import fft2, ifft2

class MultisliceWarning(UserWarning):
    pass

def ncycles(iterable, n):
    """Returns the sequence elements n times"""
    return chain.from_iterable(repeat(tuple(iterable), n))

def wdiffsim(crystal, energy, initial_wavefunction = None, resolution = (2048, 2048), 
             camera_distance = 0.2235, pixel_width = 14e-6, **kwargs):
    """
    Electron diffraction simulation in the weak-phase object approximation.

    Parameters
    ----------
    crystal : Crystal

    energy : float
        Electron energy [keV]
    initial_wavefunction : `~numpy.ndarray` or None, optional
        Initial electron wavefunction. If None (default), it will be set to a
        unit-amplitude plane-wave.
    resolution : 2-tuple of ints, optional
        Simulation resolution.
    camera_distance : float, optional
        Camera-to-sample distance [m]
    pixel_width : float, optional
        Pixel size [m].
    """
    final_resolution = resolution
    wavelength = electron_wavelength(energy)

    # Get dimensions of CCD
    vert_semiangle_max = atan( (resolution[0] * pixel_width / 2) / camera_distance)
    horz_semiangle_max = atan( (resolution[1] * pixel_width / 2) / camera_distance)
    max_k = max(vert_semiangle_max/wavelength, horz_semiangle_max/wavelength)

    X, Y, KX, KY = sim_mesh(crystal, resolution = resolution, max_k = max_k)
    if initial_wavefunction is None:
        initial_wavefunction = np.ones_like(X, dtype = np.complex)

    exit_wave = initial_wavefunction * np.exp(1j * interaction_parameter(energy) * pelectrostatic(crystal, X, Y))
    intensity = np.abs(fftshift(fft2(exit_wave)))**2
    return intensity
    # Interpolate on CCD
    #CCD_kx = np.linspace(-max_k, max_k, num = resolution[0])
    #CCD_ky = np.linspace(-max_k, max_k, num = resolution[1])
    #interpolator = RegularGridInterpolator((KX[0,:], KY[:,0]), intensity)

    #return interpolator(CCD_kx, CCD_ky)

def weak_phase(crystal, energy, initial_wavefunction = None, **kwargs):
    """
    Calculate the scattered electron wavefunction from a crystal 
    in the weak-phase object approximation

    Parameters
    ----------
    crystal : Crystal
    
    energy : float
        Electron energy [keV]
    initial_wavefunction : `~numpy.ndarray` dtype complex, or None, optional
        Initial electron wavefunction. If None (default), initial wavefunction
        is uniform plane wave.
    kwargs
        Keyword arguments are passed to simulation mesh generation
        function ``sim_mesh``
    
    Returns
    -------
    exit_wave : `~numpy.ndarray`, ndim 2, dtype complex
        Scattered electron wavefunction
    """
    X, Y, *_ = sim_mesh(crystal, **kwargs)

    if initial_wavefunction is None:
        initial_wavefunction = np.ones_like(X, dtype = np.complex)
    
    return initial_wavefunction * np.exp(1j * interaction_parameter(energy) * pelectrostatic(crystal, X, Y))

def multislice(crystal, energy, initial_wavefunction = None, diagnostic = dict(), **kwargs):
    """
    Calculate the scattered electron wavefunction from a crystal with the multislice method.

    Parameters
    ----------
    crystal : Crystal
    
    energy : float
        Electron energy [keV]
    initial_wavefunction : `~numpy.ndarray` dtype complex, or None, optional
        Initial electron wavefunction. If None (default), initial wavefunction
        is uniform plane wave.
    thickness : int, keyword-only
        Sample thickness [nm]. This parameter will be rounded to the nearest 
        unit cell dimension to avoid artifacts.
    diagnostic : dict, optional
        Dictionary that will be filled with diagnostic information.
    kwargs
        Keyword arguments are passed to simulation mesh generation
        function ``sim_mesh``
    
    Returns
    -------
    exit_wave : `~numpy.ndarray`, ndim 2, dtype complex
        Scattered electron wavefunction
    """

    X, Y, KX, KY = sim_mesh(crystal, **kwargs)

    if initial_wavefunction is None:
        initial_wavefunction = np.ones_like(X, dtype = np.complex)

    # Determine slice thickness as a divider of the unitcell period in z
    # We try to have the slice thickness be close to 2 angstroms
    # TODO: dynamically determine the perfect thickness?
    _, _, period_z = crystal.periodicity
    slice_thickness = period_z / ceil(period_z / 2)
    
    # Round thickness to the nearest unit cell
    sample_thickness = int(thickness * 10 / period_z) * period_z    # in Angs
    step_z = slice_thickness * period_z
    if step_z < 1:
        warn('A slice step of less than 1 Angstroms is not recommended', MultisliceWarning)
    
    propagator = np.exp(1j * np.pi * electron_wavelength(energy) * step_z *(KX**2 + KY**2))

    # Calculate projected electrostatic potential for each vertical slice of the unit cell
    potential_slices = list()
    min_z = min( (atom.coords[2] for atom in crystal) )
    spans = [ (step_z * i - min_z, step_z * (i + 1) - min_z) for i in range(0, int(sample_thickness/slice_thickness)) ]
    interaction = interaction_parameter(energy)
    for span_z in spans:
        integral = pelectrostatic
        potential_slices.append(np.exp(1j * interaction * pelectrostatic(crystal, X, Y, bounds = span_z)))

    wavefunction = np.array(initial_wavefunction)
    for potential_slice in ncycles(potential_slices, int(sample_thickness/period_z)):
        wavefunction = ifft2(propagator * fft2(potential_slice * wavefunction))

        # TODO: Check that relative intensity hasn't dropped too low
    
    return wavefunction

def sim_mesh(crystal, resolution = (1024, 1024), max_k = 12):
    """
    Generate real-space and frequency-space meshgrids on which to calculate
    multislice-like diffraction simulation.

    Based on the input parameters, the meshgrids' size and spacing will be adjusted
    to minimize artifacts.

    Parameters
    ----------
    crystal : Crystal instance

    resolution : 2-tuple, optional
        Desired meshgrid shape. This shape will be adjusted to minimize artifacts.
        For now, ``resolution`` must be square (e.g. ``(2048, 2048)``).
    max_k : float, optional
        Maximal simulation resolution in inverse Angstroms.
    
    Returns
    -------
    X, Y : ndarray
        Real-space meshes
    KX, KY : ndarray
        Frequency-space meshes
    """
    if resolution[0] != resolution[1]:
        raise NotImplementedError
    
    # Determine spatial frequencies on the CCD
    # From this, get real space sampling
    # Technically the spacial frequencies on the CCD are not uniformly distributed
    # but that would mess up the Fourier transforms in the Multislice method
    # Frequencies are interpolated on the CCDs later
    k_extent = np.linspace(-max_k, max_k, num = resolution[0])
    dk = k_extent[1] - k_extent[0]

    # Round real space sampling to periodicity of the lattice
    # when projected onto the x-y plane
    # WARNING: spacing in both x and y must divide a unit cell perfectly
    # Otherwise, beating frequencies will appear as artifacts.
    max_dim = np.abs(fftfreq(len(k_extent), dk)).max()
    per_x, per_y, _ = crystal.periodicity
    x_max = ceil(max_dim/per_x) * per_x
    y_max = ceil(max_dim/per_y) * per_y

    # Resolution is adjusted so that dx and dy divide per_x and per_y exactly
    # That sets the frequency range, but the frequency spacing must be also correct
    # This is done by ensuring that per_x * dx perfectly divides res_x, and same
    # for res_y
    res_x, res_y = resolution                 # initial guess
    dx, dy = 2*x_max / res_x, 2*y_max / res_y # initial guess

    dx, dy = per_x/ceil(per_x / dx), per_y/ceil(per_y / dy)

    # res = # number of unit cells * size of unit cells [px]
    res_x = ceil(res_x * dx / per_x) * (per_x / dx)
    res_y = ceil(res_y * dy / per_y) * (per_y / dy)

    # Recalculate x_max and y_max from the new resolution
    X, Y = np.meshgrid( np.arange(0, (res_x + 1) * dx, step = dx),
                        np.arange(0, (res_y + 1) * dy, step = dy))

    # Regenerate a spatial frequency sampling according to the real-space sampling
    KX, KY = np.meshgrid( fftshift(fftfreq(X.shape[1], d = X[1,1] - X[0,0])), 
                          fftshift(fftfreq(Y.shape[0], d = Y[1,1] - Y[0,0])) )

    return X, Y, KX, KY