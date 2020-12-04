# -*- coding: utf-8 -*-
"""
Robust statistics and estimators
"""

import numpy as np


def mad(arr):
    """
    Element-wise median absolute deviation (MAD) of a signal.

    .. math::

        \\text{MAD}_i = | X_i - \\text{median}(X) |

    Parameters
    ----------
    arr : array-like
        Array.

    Returns
    -------
    out : `~numpy.ndarray`, dtype float
        Median absolute difference. ``out`` will be the same shape as ``x``.
    """
    dev = np.array(arr, copy=True)
    med = np.median(dev, overwrite_input=True)
    return np.abs(arr - med)
