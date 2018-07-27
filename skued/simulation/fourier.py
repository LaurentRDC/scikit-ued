# -*- coding: utf-8 -*-
""" 
Utility functions regarding simulations and Fourier transforms.
These functions are not meant to be used outside of the simulation subpackage.
"""

import numpy as np

def fft2freq(x, y, indexing = 'xy'):
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
    if indexing == 'xy':
        extent_x, extent_y = x[0, :], y[:, 0]
    elif indexing == 'ij':
        extent_x, extent_y = x[:, 0], y[0, :]
    else:
        raise ValueError("Indexing should be either 'xy' or 'ij', not {}".format(indexing))
    
    # Spacing assuming constant x and y spacing
    spacing_x = abs(extent_x[1] - extent_x[0])
    spacing_y = abs(extent_y[1] - extent_y[0])

    freqs_x = np.fft.fftfreq(len(extent_x), d = spacing_x)
    freqs_y = np.fft.fftfreq(len(extent_y), d = spacing_y)

    return np.meshgrid(freqs_x, freqs_y, indexing = indexing)

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
    image_fft = np.fft.fft2(image)
    image_fft[K > limit] = 0.0
    return np.fft.ifft2(image_fft)