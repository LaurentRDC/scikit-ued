# -*- coding: utf-8 -*-
"""
Convenience functions for fitting time-series.
"""

import numpy as np
from math import sqrt, log
from functools import wraps
from scipy.interpolate import interp1d


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
    arr = np.full_like(time, amp + offset, dtype=float)
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
    arr = np.full_like(time, offset, dtype=float)
    arr += exponential(time, tzero, amp1, tconst1)
    arr += exponential(time, tzero, amp2, tconst2)
    return arr


def regrid(f):
    """
    Decorator that makes a function `f` evaluate correctly with uneven-spacing
    variables `t` by evaluating on a denser, even grid, and then interpolating back to `t`.

    This is useful, for example, if the function `f` involves convolutions.

    Parameters
    ----------
    f : callable
        Function of the form `f = func(t, *args, **kwargs)`, where `t` is the independent variable.

    Returns
    -------
    f_ : callable
        Callable of the form
    """
    # This function is based on an e-mail discussion with Samuel Palato. Thanks Sam!
    @wraps(f)
    def f_(time, *args, **kwargs):
        mn, mx = np.min(time), np.max(time)
        margin = mx - mn
        dt = np.abs(np.min(time[1:] - time[:-1]))
        constant_t = np.arange(mn - margin, mx + margin + dt, dt)

        y = f(constant_t, *args, **kwargs)

        intrp = interp1d(constant_t, y, kind=3, copy=False, assume_sorted=True)
        return intrp(time)

    return f_


def with_irf(fwhm):
    """
    This decorator factory that applies a Gaussian impulse response function (IRF) to a fitting function.

    Parameters
    ----------
    fwhm : float
        Full-width at half-maximum. The units of this value should be the same as the units
        used by the function it is decorating. See examples below.

    Returns
    -------
    decorator : callable
        Decorator that takes a function of the form ``f(t, *args, **kwargs)``
        and convolutes it with a Gaussian IRF.

    Examples
    --------
    Here's an example of an exponential function with an IRF. In this example,
    the ``time`` argument is in units of picoseconds. Therefore, we convolve with
    an IRF of 0.150 picoseconds (150 fs).

    >>> from skued import with_irf, exponential
    >>> @with_irf(0.150)
    ... def exponential_with_irf(time, *args, **kwargs):
    ...     return exponential(time, *args, **kwargs)

    If we were to change the definition to use femtoseconds:

    >>> from skued import with_irf, exponential # exponential is defined on picoseconds
    >>> @with_irf(150) # femtoseconds
    ... def exponential_with_irf(time, *args, **kwargs):
    ...     return exponential(time * 1000, *args, **kwargs) # defined on femtoseconds
    """

    def decorator(f):
        @wraps(f)
        @regrid
        def f_(time, *args, **kwargs):
            kernel = _gauss_kernel(time, fwhm=fwhm)
            norm = np.sum(kernel)
            return np.convolve(f(time, *args, **kwargs), kernel, mode="same") / norm

        return f_

    return decorator


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
        -((1.0 * t - t0) ** 2) / (2 * std**2)
    )
