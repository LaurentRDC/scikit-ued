# -*- coding: utf-8 -*-
"""
Multislice simulation
=====================
"""
import numpy as np
from math import ceil, fmod
from scipy.fftpack import fftfreq, fftshift

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