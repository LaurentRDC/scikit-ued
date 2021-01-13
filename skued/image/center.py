# -*- coding: utf-8 -*-
"""
Determine the center of diffraction images
==========================================
"""

import numpy as np
from skimage.registration import phase_cross_correlation
from scipy.ndimage import shift


def autocenter(im, mask=None):
    """
    Find the center of a diffraction pattern automatically.

    .. versionadded:: 2.1.1

    Parameters
    ----------
    im : ndarray, shape (N, M)
        Diffraction pattern.
    mask : ndarray, shape (N,M), dtype bool, optional
        Mask that evaluates to `True` on pixels that
        should be used to determine the center.

    Returns
    -------
    r, c : 2-tupe of ints
        Indices of the center, such that `im[r, c]` is the intensity value at
        the center of the pattern.

    Notes
    -----
    The procedure in this routine is an extension of [1]. First, the center-of-mass
    of the image intensity profile is found. This approximate center is
    then used to radially-invert the image, with the coordinate transform
    :math:`(r, \\theta) \\to (-r, \\theta)`. Diffraction patterns reflect this
    inversion symmetry, and so the shift between the original and inverted image
    is the correction to the approximate center found in the first step.

    References
    ----------
    .. [1] Liu, Lai Chung. Chemistry in Action: Making Molecular Movies with Ultrafast
           Electron Diffraction and Data Science, Chapter 2. Springer Nature, 2020.
    """

    im = np.asfarray(im)
    im -= im.min()

    if mask is None:
        mask = np.ones_like(im, dtype=np.bool)

    weights = im * mask.astype(im.dtype)

    # Center of mass. This works because the intensity envelope of a
    # diffraction pattern has a radial shape
    rr, cc = np.indices(im.shape)
    r_rough = np.average(rr, weights=weights)
    c_rough = np.average(cc, weights=weights)

    # The comparison between Friedel pairs from [1] is generalized to
    # any inversion symmetry, including polycrystalline diffraction patterns.
    im_i = radial_inversion(im, center=(r_rough, c_rough), cval=0.0)
    mask_i = radial_inversion(mask, center=(r_rough, c_rough), cval=False)

    shift = phase_cross_correlation(
        reference_image=im,
        moving_image=im_i,
        reference_mask=mask,
        moving_mask=mask_i,
    )

    return np.array([r_rough, c_rough]) + shift / 2 - np.array([1 / 2, 1 / 2])


def radial_inversion(im, center, cval):
    arr_center = np.array(im.shape) / 2
    shifted = shift(im, shift=arr_center - np.asarray(center))
    shifted = shifted[::-1, ::-1]
    return shift(
        shifted, shift=np.asarray(center) - arr_center, mode="constant", cval=cval
    )
