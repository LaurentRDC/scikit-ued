# -*- coding: utf-8 -*-
"""
Non-uniform Fast Fourier Transform.
"""

from functools import lru_cache
import numpy as np
from scipy.fftpack import fft, ifft
from math import sqrt, log, pi

def nfftfreq(M, fr = 1):
    """
    Compute the frequency range used in nfft for `M` frequency bins, resulting in the frequencies:
    
    .. math::
    
        f \in \left\{ \\text{fr} \cdot m \mid m = \overbrace{-1/2, ..., 1/2}^{M \\text{ points}} \\right\}

    Parameters
    ----------
    M : int
        Number of frequency bins.
    fr : float, optional
        Frequency range.
    
    Returns
    -------
    freqs : `~numpy.ndarray`
        Frequencies will be distributed evenly between ``-fr/2`` and ``fr/2``
    """
    M = int(M)
    return (fr/M) * np.arange(-(M // 2), M - (M // 2))

@lru_cache(maxsize = 8)
def _compute_grid_params(M, eps):
    # Choose Msp & tau from eps following Dutt & Rokhlin (1993)
    ratio = 2 if eps > 1E-11 else 3
    Msp = int(-log(eps) / (pi * (ratio - 1) / (ratio - 0.5)) + 0.5)
    Mr = max(ratio * M, 2 * Msp)
    lambda_ = Msp / (ratio * (ratio - 0.5))
    tau = pi * lambda_ / M ** 2
    return Mr, Msp, tau

def _gaussian_grid_1D(x, y, Mr, Msp, tau):
    """Compute the 1D gaussian gridding """
    ftau = np.zeros(Mr, dtype = y.dtype)
    hx = 2 * pi / Mr
    xmod = x % (2 * pi)

    m = 1 + (xmod // hx).astype(int)
    msp = np.arange(-Msp, Msp)[:, np.newaxis]
    mm = m + msp
    
    E1 = np.exp(-0.25 * (xmod - hx * m) ** 2 / tau)
    E2 = np.exp(msp * (xmod - hx * m) * pi / (Mr * tau))
    E3 = np.exp(-(pi * msp / Mr) ** 2 / tau)
    
    spread = (y * E1) * E2 * E3
    np.add.at(ftau, mm % Mr, spread)

    return ftau

def nfft(x, y, M, fr = 1.0, eps = 1E-15):
    """
    Non-uniform Fast Fourier Transform (NFFT) computed on a uniform
    frequency grid.

    Parameters
    ----------
    x : array-like
        real locations of the signal
    y : array-like
        Signal, possibly complex.
    M : int
        Number of frequencies on which the transform is computed.
    fr : float, optional
        Frequency range. Frequencies will be evenly distributed
        between ``-fr/2`` and ``fr/2``.
    eps : float, optional
        The desired approximate error for the FFT result.

    Returns
    -------
    out : `~numpy.ndarray`, dtype complex
        Non-uniform Fast Fourier Transform.

    Raises
    ------
    ValueError : if ``x`` and ``y`` don't have the same shape.

    See Also
    --------
    nfftfreq : compute the frequencies of the nfft results

    References
    ----------
    .. [NFFT] L. Greengard and J.-Y. Lee, Accelerating the Nonuniform Fast Fourier
        Transform. SIAM Review, Vol. 46, No. 3, pp. 443-454 (2005).
    """
    x, y = np.atleast_1d(x, y)
    if x.shape != y.shape:
        raise ValueError('Signal location is of unexpected shape {} \
                          compared to signal shape {} '.format(x.shape, y.shape))

    M = int(M)
    k = nfftfreq(M, fr)

    Mr, Msp, tau = _compute_grid_params(M, eps)
    ftau = _gaussian_grid_1D(x, y, Mr, Msp, tau)

    # Compute the FFT on the convolved grid
    Ftau = (1 / Mr) * fft(ftau)
    Ftau = np.concatenate([Ftau[-(M//2):], Ftau[:M//2 + M % 2]])

    # Deconvolve the grid using convolution theorem
    return (1 / len(x)) * sqrt(pi / tau) * np.exp(tau * k ** 2) * Ftau