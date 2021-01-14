# -*- coding: utf-8 -*-
"""
Determine the center of diffraction images
==========================================
"""
from os import cpu_count
import numpy as np
from skimage.registration import phase_cross_correlation
from scipy.ndimage import shift
from ..fft import with_skued_fft


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
        Indices of the center, such that ``im[r, c]`` is the intensity value at
        the center of the pattern.

    Notes
    -----
    The procedure in this routine is an extension of the one in the reference below.
    It has been adapted for both single-crystal and polycrystalline diffraction patterns
    The continuous inversion symmetry is encoded as the coordinate transformation
    :math:`(r, \\theta) \\to (-r, \\theta)`. The shift between the image and inverted image
    is the correction to the approximate center found by calculating the intensity
    center-of-mass.

    References
    ----------
    Liu, Lai Chung. Chemistry in Action: Making Molecular Movies with Ultrafast
    Electron Diffraction and Data Science, Chapter 2. Springer Nature, 2020.
    """
    if mask is None:
        mask = np.ones_like(im, dtype=np.bool)

    r_rough, c_rough = _center_of_intensity(im=im, mask=mask)

    # The comparison between Friedel pairs from [1] is generalized to
    # any inversion symmetry, including polycrystalline diffraction patterns.
    im_i = _fast_radial_inversion(im, center=(r_rough, c_rough), cval=0.0)
    mask_i = _fast_radial_inversion(mask, center=(r_rough, c_rough), cval=False)

    shift = with_skued_fft(phase_cross_correlation)(
        reference_image=im,
        moving_image=im_i,
        reference_mask=mask,
        moving_mask=mask_i,
    )

    return np.array([r_rough, c_rough]) + shift / 2 - np.array([1 / 2, 1 / 2])


def _center_of_intensity(im, mask=None):
    im = np.asfarray(im)
    im -= im.min()

    weights = im * mask.astype(im.dtype)

    rr, cc = np.indices(im.shape)
    r_rough = np.average(rr, weights=weights)
    c_rough = np.average(cc, weights=weights)
    return int(r_rough), int(c_rough)


def _fast_radial_inversion(im, center, cval):
    arr_center = np.array(im.shape) / 2
    shifted = shift(
        im, shift=arr_center - np.asarray(center), order=1, mode="constant", cval=cval
    )
    shifted = shifted[::-1, ::-1]
    return shift(
        shifted,
        shift=np.asarray(center) - arr_center,
        order=1,
        mode="constant",
        cval=cval,
    )
