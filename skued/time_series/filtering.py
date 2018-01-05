# -*- coding: utf-8 -*-

import numpy as np

def mad(arr, constant = 1.4826):
    """ 
    Median absolute deviation of a signal. 

    Parameters
    ----------
    arr : array-like
        Array. 
    constant : float, optional
        Consistency constant. Default value corresponds to normal distribution.
    
    Returns
    -------
    out : `~numpy.ndarray`, dtype float
        Median absolute difference. 
        ``out`` will be the same shape as ``x``.
    """
    arr = np.array(arr, copy = True, dtype = np.float)
    arr -= np.median(arr)
    arr[:] = np.abs(arr)
    return constant * np.median(arr)

def outliers_mad(arr, thresh = 2, constant = 1.4826):
    """
    Returns a boolean mask corresponding to outliers in an
    array, by using the median absolute deviation.

    Parameters
    ----------
    arr : array-like
        Array. 
    thresh : float, optional
        Elements of ``arr`` larger than this will be flagged
        as outliers.
    constant : float, optional
        Consistency constant. Default value corresponds to normal distribution.

    Returns
    -------
    out : `~numpy.ndarray`, dtype bool
        Outlier mask evaluates to ``True`` on elements of ``arr`` that are outliers.
    """
    out = np.zeros_like(arr, dtype = np.bool)

    med = np.median(arr)
    dev = np.abs(arr - med)
    madev = constant * np.median(dev)
    # TODO: faster to use np.logical_not(..., where = ...)?
    out[dev / madev >= thresh] = True
    return out