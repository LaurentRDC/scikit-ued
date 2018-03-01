# -*- coding: utf-8 -*-
"""
Image manipulation of powder diffraction
========================================
"""
import numpy as np
from functools import partial

from .alignment import diff_register

flip = partial(np.rot90, k = 2)

def powder_center(image, mask = None):
    """
    Finds the center of a powder diffraction pattern by comparing the
    correlation between the input and its image.

    Parameters  
    ----------
    image : `~numpy.ndarray`, ndim 2
        Image of a powder pattern
    mask : `~numpy.ndarray` or None, optional
        Mask of `image`. The mask should evaluate to `True`
        (or 1) on invalid pixels. If None (default), no mask
        is used.

    Returns
    -------
    center : 2-tuple
        Center of the powder pattern. The center is returned in the format
        relevant for array manipulations (center = [row, column] instead of 
        center = [x,y]).
    """
    if mask is None:
        mask = np.zeros_like(image, dtype = np.bool)
    mask = np.asarray(mask, dtype = np.bool)

    composite_mask = np.logical_or(mask, flip(mask))

    shift = diff_register(image, flip(image), mask = composite_mask)
    midpoints = np.array([int(axis_size / 2) for axis_size in image.shape])

    return tuple(shift[::-1]/2 + midpoints)

def _angle_bounds(bounds):
    b1, b2 = bounds
    while b1 < 0:
        b1 += 360
    while b1 > 360:
        b1 -= 360
    while b2 < 0:
        b2 += 360
    while b2 > 360:
        b2 -= 360
    return tuple(sorted((b1, b2)))
    
def azimuthal_average(image, center, mask = None, angular_bounds = None, trim = True):
    """
    This function returns an azimuthally-averaged pattern computed from an image, 
    e.g. polycrystalline diffraction.

    Parameters
    ----------
    image : array_like, shape (M, N)
        Array or image.
    center : array_like, shape (2,)
        coordinates of the center (in pixels).
    mask : `~numpy.ndarray` or None, optional
        Evaluates to True on invalid elements of array.
    angular_bounds : 2-tuple or None, optional
        If not None, the angles between first and second elements of `angular_bounds`
        (inclusively) will be used for the average. Angle bounds are specified in degrees.
        0 degrees is defined as the positive x-axis. Angle bounds outside [0, 360) are mapped back
        to [0, 360).
    trim : bool, optional
        If True, leading and trailing zeros (possible due to the usage of masks) are trimmed.

    Returns
    -------
    radius : `~numpy.ndarray`, ndim 1
        Radius of the average [px]. ``radius`` might not start at zero, depending on the ``trim`` parameter.
    average : `~numpy.ndarray`, ndim 1
        Angular-average of the array.
    """
    if mask is None:
        mask = np.zeros_like(image, dtype = np.bool)

    xc, yc = center  

    # Create meshgrid and compute radial positions of the data
    # The radial positions are rounded to the nearest integer
    # TODO: interpolation? or is that too slow?
    Y, X = np.indices(image.shape)
    R = np.hypot(X - xc, Y - yc)
    Rint = np.rint(R).astype(np.int)

    if angular_bounds:
        mi, ma = _angle_bounds(angular_bounds)
        angles = np.rad2deg(np.arctan2(Y - yc, X - xc)) + 180  # arctan2 is defined on [-pi, pi] but we want [0, pi]
        in_bounds = np.logical_and(mi <= angles, angles <= ma)
    else:
        in_bounds = np.ones_like(image, dtype = np.bool)

    valid = np.logical_not(mask)[in_bounds]
    image = image[in_bounds]
    Rint = Rint[in_bounds]

    px_bin = np.bincount(Rint, weights = valid*image)
    r_bin = np.bincount(Rint, weights = valid)
    radius = np.arange(0, r_bin.size)

    # Make sure r_bin is never 0 since it it used for division anyway
    np.maximum(r_bin, 1, out = r_bin)

    # We ignore the leading and trailing zeroes, which could be due to
    first, last = 0, -1
    if trim:
        first, last = _trim_bounds(px_bin)

    return radius[first:last], px_bin[first:last]/r_bin[first:last]

def _trim_bounds(arr):
    """ Returns the bounds which would be used in numpy.trim_zeros """
    first = 0
    for i in arr:
        if i != 0.:
            break
        else:
            first = first + 1
    last = len(arr)
    for i in arr[::-1]:
        if i != 0.:
            break
        else:
            last = last - 1
    return first, last