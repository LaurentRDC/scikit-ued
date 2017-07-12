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
def nfold(im, mod, center = None, mask = None, **kwargs):
    """ 
    Returns an images averaged according to n-fold rotational symmetry.
    Keyword arguments are passed to skimage.transform.rotate()

    Parameters
    ----------
    im : array_like, ndim 2
        Image to be averaged.
    center : array_like, shape (2,) or None, optional
        coordinates of the center (in pixels). If ``center=None``, the image is rotated around
        its center, i.e. ``center=(rows / 2 - 0.5, cols / 2 - 0.5)``.
    mod : int
        Fold symmetry number. Valid numbers must be a divisor of 360.
    mask : `~numpy.ndarray` or None, optional
        Mask of `image`. The mask should evaluate to `True`
        (or 1) on invalid pixels. If None (default), no mask
        is used.
    

    Returns
    -------
    out : `~numpy.ndarray`, dtype float
        Averaged image.

    Raises
    ------
    ValueError
        If `mod` is not a divisor of 360 deg.
    
    See also
    --------
    skimage.transform.rotate : Rotate images by interpolation.
    """
    if (360 % mod) != 0:
        raise ValueError('{}-fold rotational symmetry is not valid (not a divisor of 360).'.format(mod))
    angles = range(0, 360, int(360/mod))

    im = np.array(im, dtype = np.float, copy = True)

    if mask is not None:
       im[mask] = np.nan
    
    rotate_kwargs = {'mode': 'constant', 'preserve_range': True}
    rotate_kwargs.update(kwargs)

    stack = np.dstack([rotate(im, angle, center = center, **rotate_kwargs) for angle in angles])
    avg = np.nanmean(stack, axis = 2)
    return np.nan_to_num(avg)

def nfold_symmetry(*args, **kwargs):
    warn('nfold_symmetry() is deprecated. Please use nfold in \
          the future, as it supports more features.', DeprecationWarning)
    return nfold(*args, **kwargs)
