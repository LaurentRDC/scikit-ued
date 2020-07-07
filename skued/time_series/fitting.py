# -*- coding: utf-8 -*-
"""
Convenience functions for fitting time-series.
"""

import numpy as np


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
