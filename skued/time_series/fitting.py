# -*- coding: utf-8 -*-
"""
Convenience functions for fitting time-series.
"""

import numpy as np

def exponential_decay(time, tzero, amp, tconst, offset = 0):
    """
    Exponential decay curve with onset. Equivalent to the following function:

    .. math::

        I(t) =
        \\begin{cases} 
            I_0 + O &\\text{if } t < t_0 \\\\
            I_0 e^{-\lambda (t - t_0)} + O &\\text{if } t \ge t_0
        \\end{cases}

    Parameters
    ----------
    time : `~numpy.ndarray`
        Time values array [ps].
    tzero : float
        Time-zero :math:`t_0`. Exponential decay happens for `time > tzero`.
    amp : float
        Initial amplitude :math:`I_0`.
    tconst : float
        Decay time-constant :math:`\lambda`.
    offset : float, optional
        Constant amplitude offset :math:`O`.
    """
    arr = np.full_like(time, amp + offset, dtype = np.float)
    arr[time > tzero] = amp * np.exp(-tconst * (time[time > tzero] - tzero)) + offset
    return arr
