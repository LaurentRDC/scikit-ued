# -*- coding: utf-8 -*-
"""
Convenience functions for fitting time-series.
"""

import numpy as np

# TODO: biexponential
#       triexponential
#       oscillation after onset

def exponential_decay(time, tzero, amp, tconst, offset = 0):
    """
    Exponential decay curve with onset. Equivalent to the following function:

    .. math::

        I(t) =
        \\begin{cases} 
            I_0 + O &\\text{if } t < t_0 \\\\
            I_0 e^{-(t - t_0)/\\tau} + O &\\text{if } t \ge t_0
        \\end{cases}

    This functions is expected to be used in conjunction with 
    ``scipy.optimize.curve_fit`` or similar fitting routines.

    Parameters
    ----------
    time : `~numpy.ndarray`, shape(N,)
        Time values array [ps].
    tzero : float
        Time-zero :math:`t_0` [ps]. Exponential decay happens for :math:`t > t_0`.
    amp : float
        Initial amplitude :math:`I_0`.
    tconst : float
        Decay time-constant :math:`\\tau` [ps].
    offset : float, optional
        Constant amplitude offset :math:`O`.

    Returns
    -------
    exp_decay : `~numpy.ndarray`, shape (N,)
        Exponential decay curve.
    
    See also
    --------
    scipy.optimize.curve_fit : 1D curve-fitting routine
    """
    arr = np.full_like(time, amp + offset, dtype = np.float)
    arr[time > tzero] = amp * np.exp(-(time[time > tzero] - tzero)/tconst) + offset
    return arr
