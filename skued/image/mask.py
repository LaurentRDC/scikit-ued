# -*- coding: utf-8 -*-
"""
Image mask routines
===================
"""
from collections import Iterable
from itertools import tee

import numpy as np

from npstreams import array_stream, istd, peek, last, iprod

# array_stream decorator ensures that input images are cast to ndarrays
@array_stream
def mask_from_collection(images, px_thresh = (0, 3e4), std_thresh = None):
    """ 
    Determine binary mask from a set of images. These images should represent identical measurements, e.g. a set
    of diffraction patterns before photoexcitation. Pixels are rejected on the following two criteria:

        * Pixels with a value above a certain threshold or below zero, for any image in the set, are considered dead;
        * Pixels with a cumulative standard deviation above a certain threshold are considered uncertain.

    This function operates in constant memory, so large collections of images can be used.

    Parameters
    ----------
    images : iterable of ndarray
        These images should represent identical measurements. ``images`` can also be a generator.
    px_thresh : float or iterable, optional
        Pixels with a value outside of [``min(px_thresh)``, ``max(px_thresh)``] in any of the images in ``images``
        are rejected. If ``px_thresh`` is a single float, it is assumed to be the maximal intensity, and no lower
        bound is enforced.
    std_thresh : int or float or None, optional
        Standard-deviation threshold. If the standard deviation of a pixel exceeds ``std_thresh``, 
        it is rejected. If None (default), a threshold is not enforced.
    
    Returns
    -------
    mask : `~numpy.ndarray`, dtype bool
        Pixel mask. Pixels where ``mask`` is True are invalid.
    
    Notes
    -----
    ``numpy.inf`` can be used to have a lower pixel value bound but no upper bound. For example, to
    reject all negative pixels only, set ``px_thresh = (0, numpy.inf)``
    """
    if isinstance(px_thresh, Iterable):
        min_int, max_int = min(px_thresh), max(px_thresh)
    else:
        min_int, max_int = None, px_thresh
    
    first, images = peek(images)
    mask = np.zeros_like(first, dtype = np.bool)    # 0 = False

    images, images_for_std = tee(images)
    for image, std in zip(images, istd(images_for_std)):
        
        mask[image > max_int] = True
        
        if std_thresh is not None:
            mask[std > std_thresh] = True

        if min_int is not None:
            mask[image < min_int] = True
    
    return mask

def combine_masks(*masks):
    """ 
    Combine multiple pixel masks into one. This assumes that pixel masks evaluate
    to ``True`` on invalid pixels.

    Returns
    -------
    combined : `~numpy.ndarray`, dtype bool 
    """
    # By multiplying boolean arrays, values of False propagate
    # Hence, much easier to do if invalue pixels are False instead of True
    valids = map(np.logical_not, masks)
    combined_valid = last(iprod(valids, dtype = np.bool))
    return np.logical_not(combined_valid)

def mask_image(image, mask, fill_value = 0, copy = True):
    """
    Fill invalid pixels in an image with another value, according to a pixel mask.
    While this function has simply functionality, it is ideal for integrating into a
    pipeline with ``functools.partial``.

    Parameters
    ---------
    image : `~numpy.ndarray`
    
    mask : `~numpy.ndarray`
        Boolean array. ``mask`` should evaluate to ``True`` on invalid pixels.
    fill_value : float, optional
        Invalid pixels fill value.
    copy : bool, optional
        If True (default), ``image`` is copied before masking. If False, ``image`` is modified in-place.
    
    Returns
    -------
    masked : `~numpy.ndarray`
        Masked image. If ``copy = True``, masked points to the same object as ``image``.
    """
    image = np.array(image, copy = copy)
    mask = np.asarray(mask, dtype = np.bool)

    image[mask] = fill_value
    return image