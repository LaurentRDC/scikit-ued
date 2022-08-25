# -*- coding: utf-8 -*-
"""
Non-uniform Fast Fourier Transform.
"""

from math import log, pi, sqrt

import numpy as np
from scipy.fftpack import fft


def nfftfreq(M, df=1):
    """
    Compute the frequency range used in nfft for `M` frequency bins.

    Parameters
    ----------
    M : int
        Number of frequency bins.
    df : float, optional
        Frequency range.

    Returns
    -------
    freqs : `~numpy.ndarray`
    """
    M = int(M)
    return df * np.arange(-(M // 2), M - (M // 2))


def nfft(x, y, M, df=1.0, eps=1e-15):
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
    df : float, optional
        Frequency range.
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
        raise ValueError(
            f"Signal location is of unexpected shape {x.shape} \
                          compared to signal shape {y.shape} "
        )

    M = int(M)
    k = nfftfreq(M, df)

    R = 3
    Mr = R * M
    Msp = int(-log(eps) / (pi * (R - 1) / (R - 0.5)) + 0.5)
    tau = pi * Msp / (R * (R - 0.5)) / M**2
    ftau = _fast_gaussian_grid(x, y, Mr, Msp, tau)

    # Compute the FFT on the convolved grid
    Ftau = (1 / Mr) * fft(ftau)
    Ftau = np.concatenate([Ftau[-(M // 2) :], Ftau[: M // 2 + M % 2]])

    # Deconvolve the grid using convolution theorem
    return (1 / len(x)) * sqrt(pi / tau) * np.exp(tau * k**2) * Ftau


def _fast_gaussian_grid(x, y, Mr, Msp, tau):
    ftau = np.zeros(Mr, dtype=y.dtype)
    hx = 2 * pi / Mr
    xmod = x % (2 * pi)

    m = 1 + (xmod // hx).astype(int)
    msp = np.arange(-Msp, Msp)[:, np.newaxis]
    mm = m + msp

    E1 = np.exp(-0.25 * (xmod - hx * m) ** 2 / tau)
    E2 = np.exp(msp * (xmod - hx * m) * pi / (Mr * tau))
    E3 = np.exp(-((pi * msp / Mr) ** 2) / tau)

    spread = (y * E1) * E2 * E3
    np.add.at(ftau, mm % Mr, spread)

    return ftau
