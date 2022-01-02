# -*- coding: utf-8 -*-
"""
Kinematic simulation of diffraction patterns
============================================
"""

import numpy as np
import scipy.fft as fft
from ..fft import with_skued_fft
from .potential import pelectrostatic
from ..eproperties import interaction_parameter
from scipy.interpolate import RegularGridInterpolator


@with_skued_fft
def kinematicsim(crystal, kx, ky, energy=90):
    """
    Propagate a plane wave through a crystal and compute the resulting
    diffraction pattern, in the kinematic approximation (thin specimen).

    .. versionadded:: 2.0.5

    Parameters
    ----------
    crystal : crystals.Crystal
        Crystal from which to scatter.
    kx, ky :  `~numpy.ndarray`, shape (N,M)
        Momenta mesh where to calculate the diffraction pattern [:math:`Å^{-1}`]
    energy : float, optional
        Electron energy [keV]

    Returns
    -------
    diff_pattern : `~numpy.ndarray`
        Scattered intensity.
    """
    shape = tuple(map(fft.next_fast_len, kx.shape))
    period_x, period_y, period_z = crystal.periodicity

    # We create the grid ourselves so that we minimize Fourier artifacts as much as possible.
    # It is much easier to interpolate to the requested grid than to prevent artifact formation.
    extent = 8 * period_x * period_y
    extent_x = np.linspace(0, extent, num=shape[0])
    extent_y = np.linspace(0, extent, num=shape[1])

    xx, yy = np.meshgrid(
        extent_x,
        extent_y,
        indexing="xy",
    )
    kx_, ky_ = fft2freq(xx, yy, indexing="xy")
    k = np.hypot(kx_, ky_)

    potential = pelectrostatic(crystal, xx, yy)
    transmission_function = np.exp(1j * interaction_parameter(energy) * potential)

    exit_wave = fft.ifft2(
        fft.fft2(np.ones_like(xx, dtype=complex) * transmission_function)
    )
    intensity = fft.fftshift(np.abs(fft.fft2(exit_wave)) ** 2)

    kx_ = fft.fftshift(kx_)
    ky_ = fft.fftshift(ky_)

    # Note that the definition of 'frequency' in fftfreq & friends necessitates dividing by 2pi
    twopi = 2 * np.pi
    return RegularGridInterpolator(
        points=(kx_[0, :], ky_[:, 0]),
        values=intensity,
        bounds_error=False,
        fill_value=0,
    ).__call__(xi=(kx / twopi, ky / twopi))


def fft2freq(x, y, indexing="xy"):
    """
    Return the Discrete Fourier Transform sample frequencies for a 2D array defined on ``x`` and ``y``.
    Generalization of ``fftfreq``.

    Parameters
    ----------
    x, y : `~numpy.ndarray`, ndim 2
        Meshgrid-style arrays. Spacing must be uniform.
    indexing : {'ij', 'xy'}, optional
        Indexing used to generate ``x`` and ``y``.

    Returns
    -------
    kx, ky : `~numpy.ndarray`, ndim 2

    Raises
    ------
    ValueError : if ``indexing`` is invalid.
    """
    if indexing == "xy":
        extent_x, extent_y = x[0, :], y[:, 0]
    elif indexing == "ij":
        extent_x, extent_y = x[:, 0], y[0, :]
    else:
        raise ValueError(
            "Indexing should be either 'xy' or 'ij', not {}".format(indexing)
        )

    # Spacing assuming constant x and y spacing
    spacing_x = abs(extent_x[1] - extent_x[0])
    spacing_y = abs(extent_y[1] - extent_y[0])

    freqs_x = fft.fftfreq(len(extent_x), d=spacing_x)
    freqs_y = fft.fftfreq(len(extent_y), d=spacing_y)

    return np.meshgrid(freqs_x, freqs_y, indexing=indexing)


def limit_bandwidth(image, K, limit):
    """
    Limit the bandwidth of an image.

    Parameters
    ----------
    image : `~numpy.ndarray`, ndim 2
        Image to be bandwidth-limited.
    K : `~numpy.ndarray`, ndim 2
        Wavevector norm on which the Fourier transform of the image is defined.
    limit : float
        Bandwidth limit.

    Returns
    -------
    limited : `~numpy.ndarray`
        Bandwidth-limited image.
    """
    image_fft = fft.fft2(image)
    image_fft[K > limit] = 0.0
    return fft.ifft2(image_fft)
