# -*- coding: utf-8 -*-
"""
Image manipulation involving symmetry
"""
from functools import partial
from skimage.transform import rotate
import numpy as np

# TODO: out parameter?
def nfold_symmetry(im, center, mod, **kwargs):
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

    im = np.asarray(im)
    angles = range(0, 360, int(360/mod))

    kwargs.update({'preserve_range': True})
    return sum(rotate(im, angle, center = center, **kwargs) for angle in angles)/len(angles)