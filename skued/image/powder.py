# -*- coding: utf-8 -*-
"""
Image manipulation of powder diffraction
========================================
"""
import numpy as np
from warnings import warn
from .alignment import diff_register

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

	shift = diff_register(image, np.rot90(image, k = 2), mask = mask * np.rot90(mask, k = 2))
	midpoints = np.array([int(axis_size / 2) for axis_size in image.shape])
	
	# I have found that there is always a residual (0.5, 0.5)
	center = shift[::-1]/2 + midpoints - np.array([0.5, 0.5])
	return tuple(center)

def angular_average(*args, **kwargs):
    warn('angular_average is deprecated. See azimuthal_average for more features.', DeprecationWarning)
    return azimuthal_average(*args, **kwargs)

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
    
def azimuthal_average(image, center, mask = None, angular_bounds = None):
    """
    This function returns an angularly-averaged pattern computed from a diffraction pattern, 
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

    Returns
    -------
    radius : `~numpy.ndarray`, ndim 1
        Radius of the average [px]. If ``mask`` is provided, ``radius`` might not start at one;
        ``average`` is trimmed of leading invalid pixels.
    average : `~numpy.ndarray`, ndim 1
        Angular-average of the array.
    """
    # TODO: interpolation?
    # TODO: error?
    if mask is None:
        mask = np.zeros_like(image, dtype = np.bool)

    xc, yc = center  

    #Create meshgrid and compute radial positions of the data
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

    # We ignore the leading and trailing zeroes, which could be due to
    first, last = _trim_bounds(px_bin)
    radial_intensity = px_bin[first:last]/r_bin[first:last]

    # Error as the standard error in the mean, at each pixel
    # Standard error = std / sqrt(N)
    # std = sqrt(var - mean**2)
    #if extras is not None:
    #    var_bin = np.bincount(R, weights = image**2)[first:last]/r_bin[first:last]
    #    radial_intensity_error = np.sqrt(var_bin - radial_intensity**2)/np.sqrt(r_bin[first:last])

    return radius[first:last], radial_intensity

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