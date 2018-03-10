# -*- coding: utf-8 -*-
"""
Image manipulation involving symmetry
=====================================
"""
import numpy as np
from skimage.transform import rotate

from npstreams import average, nan_to_num


# TODO: out parameter?
def nfold(im, mod, center = None, mask = None, fill_value = 0.0):
    """ 
    Returns an images averaged according to n-fold rotational symmetry. This can be used to
    boost the signal-to-noise ratio on an image with known symmetry, e.g. a diffraction pattern.

    Parameters
    ----------
    im : array_like, ndim 2
        Image to be azimuthally-symmetrized.
    center : array_like, shape (2,) or None, optional
        Coordinates of the center (in pixels). If ``center=None``, the image is rotated around
        its center, i.e. ``center=(rows / 2 - 0.5, cols / 2 - 0.5)``.
    mod : int
        Fold symmetry number. Valid numbers must be a divisor of 360.
    mask : `~numpy.ndarray` or None, optional
        Mask of `image`. The mask should evaluate to `True` (or 1) on invalid pixels. 
        If None (default), no mask is used.
    fill_value : float, optional
        In the case of a mask that overlaps with itself when rotationally averaged,
        the overlapping regions will be filled with this value.

    Returns
    -------
    out : `~numpy.ndarray`, dtype float
        Symmetrized image.

    Raises
    ------
    ValueError : If `mod` is not a divisor of 360 deg.
    """
    if (360 % mod):
        raise ValueError('{}-fold rotational symmetry is not valid (not a divisor of 360).'.format(mod))
    angles = range(0, 360, int(360/mod))

    # Data-type must be float because of use of NaN
    im = np.array(im, dtype = np.float, copy = True)

    if mask is not None:
       im[mask] = np.nan
    
    kwargs = {'center': center, 'mode' : 'constant', 'cval' : 0, 'preserve_range': True}

    # Use weights because edges of the pictures, which might be cropped by the rotation
    # should not count in the average
    wt = np.ones_like(im, dtype = np.uint8)
    weights = (rotate(wt, angle, **kwargs) for angle in angles)
    rotated = (rotate(im, angle, **kwargs) for angle in angles)

    avg = average(rotated, weights = weights, ignore_nan = True)
    return nan_to_num(avg, fill_value, copy = False)