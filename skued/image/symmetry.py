# -*- coding: utf-8 -*-
"""
Image manipulation involving symmetry
=====================================
"""
from functools import partial, wraps
from skimage.transform import rotate
import numpy as np
from warnings import warn

# TODO: out parameter?
def nfold(im, center, mod, mask = None, **kwargs):
    """ 
    Returns an images averaged according to n-fold rotational symmetry.
    Keyword arguments are passed to skimage.transform.rotate()

    Parameters
    ----------
    im : array_like, ndim 2
        Image to be averaged.
    center : array_like, shape (2,)
        coordinates of the center (in pixels).
    mod : int
        Fold symmetry number. Valid numbers must be a divisor of 360.
    mask : `~numpy.ndarray` or None, optional
        Mask of `image`. The mask should evaluate to `True`
        (or 1) on invalid pixels. If None (default), no mask
        is used.

    Returns
    -------
    out : `~numpy.ndarray`

    Raises
    ------
    ValueError
        If `mod` is not a divisor of 360 deg.
    """
    if (360 % mod) != 0:
        raise ValueError('Rotational symmetry of {} is not valid.'.format(mod))

    if mask is None:
        mask = np.zeros_like(im, dtype = np.bool)

    im = np.array(im, dtype = np.float, copy = True)
    im[mask] = np.nan
    angles = range(0, 360, int(360/mod))

    kwargs.update({'preserve_range': True})
    stack = np.dstack([rotate(im, angle, center = center, **kwargs) for angle in angles])

    avg = np.nanmean(stack, axis = 2)
    return np.nan_to_num(avg)

def nfold_symmetry(*args, **kwargs):
    warn('nfold_symmetry() is deprecated. Please use nfold in \
          the future, as it supports more features.', DeprecationWarning)
    return nfold(*args, **kwargs)
