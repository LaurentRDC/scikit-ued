# -*- coding: utf-8 -*-
"""
Convenience functions for fitting time-series.
"""

import numpy as np
from math import sqrt, log


def exponential(time, tzero, amp, tconst, offset=0):
    """
    Exponential curve with onset. Equivalent to the following function:

    .. math::

        I(t) = I_0 (1 - e^{-(t - t_0)/\\tau}) \\Theta(t - t_0) + O

    This functions is expected to be used in conjunction with 
    ``scipy.optimize.curve_fit`` or similar fitting routines.

    Parameters
    ----------
    time : `~numpy.ndarray`, shape(N,)
        Time values array [ps].
    tzero : float
        Time-zero :math:`t_0` [ps]. Exponential behavior happens for :math:`t > t_0`.
    amp : float
        Initial amplitude :math:`I_0`.
    tconst : float
        Time-constant :math:`\\tau` [ps].
    offset : float, optional
        Constant amplitude offset :math:`O`.

    Returns
    -------
    exp : `~numpy.ndarray`, shape (N,)
        Exponential curve.
    
    See also
    --------
    biexponential : bi-exponential curve with onset
    """
    return (
        np.heaviside(time - tzero, 1 / 2) * (amp * (1 - np.exp(-(time - tzero) / tconst)))
        + offset
    )


def biexponential(time, tzero, amp1, amp2, tconst1, tconst2, offset=0):
    """
    Bi-exponential curve with onset. Equivalent to the following function:

    .. math::

        I(t) = \\Theta(t - t_0) \\left[ I_1 (1 - e^{-(t - t_0)/\\tau_1}) + I_2 (1 - e^{-(t - t_0)/\\tau_2})\\right] + O

    This functions is expected to be used in conjunction with 
    ``scipy.optimize.curve_fit`` or similar fitting routines.

    Parameters
    ----------
    time : `~numpy.ndarray`, shape(N,)
        Time values array [ps].
    tzero : float
        Time-zero :math:`t_0` [ps]. Exponential behavior happens for :math:`t > t_0`.
    amp1 : float
        Initial amplitude :math:`I_1`.
    amp2 : float
        Initial amplitude :math:`I_2`.
    tconst1 : float
        Decay time-constant :math:`\\tau_1` [ps].
    tconst2 : float
        Decay time-constant :math:`\\tau_2` [ps].
    offset : float, optional
        Constant amplitude offset :math:`O`.

    Returns
    -------
    biexp : `~numpy.ndarray`, shape (N,)
        Biexponential curve.
    
    See also
    --------
    exponential : single-exponential curve with onset
    """
    arr = np.full_like(time, offset, dtype=np.float)
    arr += exponential(time, tzero=tzero, amp=amp1, tconst=tconst1, offset=0)
    arr += exponential(time, tzero=tzero, amp=amp2, tconst=tconst2, offset=0)
    return arr


def gauss_kernel(t, t0, fwhm):
    """
    Gaussian convolution kernel.

    Parameters
    ----------
    x : array-like
        Independant variable
    t0 : array-like
        t0 offset
    fwhm : float
        Full-width at half-maximum of the Gaussian kernel
    """
    std = fwhm / (2 * sqrt(2 * log(2)))
    return (1/(np.sqrt(2*np.pi)*std)) * np.exp(-(1.0*t-t0)**2 / (2*std**2))


def convolve(arr, kernel):
    """ Convolution of array with kernel. """
    #logger.debug("Convolving...")
    npts = min(len(arr), len(kernel))
    pad  = np.ones(npts)
    tmp  = np.concatenate((pad*arr[0], arr, pad*arr[-1]))
    norm = np.sum(kernel)
    out  = np.convolve(tmp, kernel, mode='same')
    noff = int((len(out) - npts)/2)
    return out[noff:noff+npts]/norm

# TODO: test with unevenly-spaced data points
# TODO: figure out wth is going on with the width
def with_irf(fwhm, f):
    """
    This decorator applies a Gaussian impulse response function (IRF) to a fitting function.

    Parameters
    ----------
    fwhm : float
        Full-width at half-maximum.
    f : callable
        Fit function (e.g. :func:`exponential`)
    
    Returns
    -------
    f_ : callable
        Transformed function
    
    Examples
    --------

    """
    def f_(time, *args, **kwargs):
        kernel = gauss_kernel(time, t0=0, fwhm=fwhm)
        return convolve(f(time, *args, **kwargs), kernel)
    return f_