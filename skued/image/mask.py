# -*- coding: utf-8 -*-
"""
Image mask routines
===================
"""
from collections import Iterable
from itertools import tee
import numpy as np
from npstreams import peek, array_stream, istd

# array_stream decorator ensures that input images are cast to ndarrays
@array_stream
def mask_from_collection(images, px_thresh = (0, 3e4), std_thresh = 10):
    """ 
    Determine binary mask from a set of images. These images should represent identical measurements, e.g. a set
    of diffraction patterns before photoexcitation. Pixels are rejected on the following two criteria:

        * Pixels with a value above a certain threshold or below zero, for any image in the set, are considered dead;
        * Pixels with a cumulative standard deviation above a certain threshold are considered uncertain.

    Parameters
    ----------
    images : iterable of ndarray
        These images should represent identical measurements. ``images`` can also be a generator.
    px_thresh : float or iterable, optional
        Pixels with a value outside of [``min(px_thresh)``, ``max(px_thresh)``] in any of the images in ``images``
        are rejected. If ``px_thresh`` is a single float, it is assumed to be the maximal intensity, and no lower
        bound is enforced.
    std_thresh : int or float, optional
        Standard-deviation threshold. If the standard deviation of a pixel exceeds ``std_thresh``, 
        it is rejected.
    
    Returns
    -------
    mask : `~numpy.ndarray`, dtype bool
        Pixel mask. Pixels where ``mask`` is True are invalid.
    """
    if isinstance(px_thresh, Iterable):
        min_int = min(px_thresh)
        max_int = max(px_thresh)
    else:
        min_int = None
        max_int = px_thresh
    
    first, images = peek(images)
    mask = np.zeros_like(first, dtype = np.bool)    # 0 = False

    images, images_for_std = tee(images)
    for image, std in zip(images, istd(images_for_std)):
        
        mask[image > max_int] = True
        mask[std > std_thresh] = True

        if min_int is not None:
            mask[image < min_int] = True
    
    return mask