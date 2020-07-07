# -*- coding: utf-8 -*-
"""
Convenience functions for fitting time-series.
"""

import numpy as np
from math import sqrt, log
from functools import wraps


def exponential(time, tzero, amp, tconst, offset=0):
    """
    Exponential curve with onset. Equivalent to the following function:

    .. math::

        I(t) =
        \\begin{cases} 
            I_0 + O &\\text{if } t < t_0 \\\\
            I_0 e^{-(t - t_0)/\\tau} + O &\\text{if } t \ge t_0
        \\end{cases}

    where :math:`\\Theta(t - t_0)` is the Heaviside step function. 

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
        Decay time-constant :math:`\\tau` [ps].
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
    arr = np.full_like(time, amp + offset, dtype=np.float)
    arr[time > tzero] = amp * np.exp(-(time[time > tzero] - tzero) / tconst) + offset
    return arr


def biexponential(time, tzero, amp1, amp2, tconst1, tconst2, offset=0):
    """
    Bi-exponential curve with onset. Equivalent to the following function:

    .. math::

        I(t) =
        \\begin{cases} 
            I_1 + I_2 + O &\\text{if } t < t_0 \\\\
            I_1 e^{-(t - t_0)/\\tau_1} + I_2 e^{-(t - t_0)/\\tau_2} + O &\\text{if } t \ge t_0
        \\end{cases}

    where :math:`\\Theta(t - t_0)` is the Heaviside step function.

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
    arr += exponential(time, tzero, amp1, tconst1)
    arr += exponential(time, tzero, amp2, tconst2)
    return arr


# TODO: test with unevenly-spaced data points
# TODO: add example
# TODO: add example to user guide
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

    @wraps(f)
    def f_(time, *args, **kwargs):
        kernel = _gauss_kernel(time, fwhm=fwhm)
        return _convolve(f(time, *args, **kwargs), kernel)

    return f_


def _gauss_kernel(t, fwhm):
    """
    Gaussian convolution kernel.

    Parameters
    ----------
    t : array-like
        Independant variable
    fwhm : float
        Full-width at half-maximum of the Gaussian kernel
    """
    # It is important that the kernel be centered in the array `t`
    # for the convolution operation to work properly
    t0 = t[int(len(t) / 2)]

    std = fwhm / (2 * sqrt(2 * log(2)))
    return (1 / (np.sqrt(2 * np.pi) * std)) * np.exp(
        -((1.0 * t - t0) ** 2) / (2 * std ** 2)
    )


def _convolve(arr, kernel):
    """ Convolution of array with kernel. """
    npts = min(len(arr), len(kernel))
    pad = np.ones(npts)
    tmp = np.concatenate((pad * arr[0], arr, pad * arr[-1]))
    norm = np.sum(kernel)
    out = np.convolve(tmp, kernel, mode="same")
    noff = int((len(out) - npts) / 2)
    return out[noff : noff + npts] / norm
