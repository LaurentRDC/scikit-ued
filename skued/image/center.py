# -*- coding: utf-8 -*-
"""
Determine the center of diffraction images
==========================================
"""

from math import floor
import numpy as np
from skimage.registration import phase_cross_correlation
from scipy.ndimage import gaussian_filter
from ..fft import with_skued_fft
from warnings import catch_warnings, simplefilter


def autocenter(im, mask=None, normalize_bg=True):
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
    normalize_bg: bool, optional
        If `True` (default), an attempt will be made to remove
        asymmetries in the background of `im`. This can sometimes
        provide worse results, so you may want to disable this feature.

        .. versionadded:: 2.1.16


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
    determines the correction to the approximate center found by calculating the intensity
    center-of-mass.

    References
    ----------
    Liu, Lai Chung. Chemistry in Action: Making Molecular Movies with Ultrafast
    Electron Diffraction and Data Science, Chapter 2. Springer Nature, 2020.
    """

    im = np.array(im, copy=True, dtype=float)
    im -= im.min()

    if mask is None:
        mask = np.ones_like(im, dtype=bool)
        weights = im
    else:
        weights = im * mask.astype(im.dtype)

    rr, cc = np.indices(im.shape)
    r_ = int(np.average(rr, weights=weights))
    c_ = int(np.average(cc, weights=weights))

    # Determine the smallest center -> side distance, and crop around that
    # This is for two reasons.
    # 1. Some diffraction patterns are not centered, and so there's a lot
    # of image area that cannot be used for registration.
    # 2. radial inversion becomes simple inversion of dimensions
    side_length = floor(min([r_, abs(r_ - im.shape[0]), c_, abs(c_ - im.shape[1])]))
    rs = slice(r_ - side_length, r_ + side_length)
    cs = slice(c_ - side_length, c_ + side_length)
    im = im[rs, cs]
    mask = mask[rs, cs]

    # Certain images display a gradient in the overall intensity of diffraction
    # peaks that come from ewald sphere walkoff
    # e.g. (n00) systematically brighter than (-n00)
    # For this purpose, we normalize the intensity by some "background",
    # i.e. very blurred diffraction pattern
    # This step is optional because it may negatively affect the results
    # on very clean images. See issue #45
    # (https://github.com/LaurentRDC/scikit-ued/issues/45#issuecomment-2180808898)
    if normalize_bg:
        with catch_warnings():
            simplefilter("ignore", category=RuntimeWarning)
            im /= gaussian_filter(input=im, sigma=min(im.shape) / 25, truncate=2)
    im = np.nan_to_num(im, copy=False)

    # The comparison between Friedel pairs from [1] is generalized to
    # any inversion symmetry, including polycrystalline diffraction patterns.
    im_i = im[::-1, ::-1]
    mask_i = mask[::-1, ::-1]

    # masked normalized cross-correlation is extremely expensive
    # we therefore downsample large images for essentially identical result
    # but ~4x decrease in processing time
    downsampling = 1
    if min(im.shape) > 1024:
        downsampling = 2

    shift, *_ = with_skued_fft(phase_cross_correlation)(
        reference_image=im[::downsampling, ::downsampling],
        moving_image=im_i[::downsampling, ::downsampling],
        reference_mask=mask[::downsampling, ::downsampling],
        moving_mask=mask_i[::downsampling, ::downsampling],
    )
    # Because images were downsampled, the correction
    # factor to the rough center should be increased from the measured shift
    correction = shift * downsampling

    return np.array([r_, c_]) + correction / 2


def auto_masking(im, threshold=0.1):
    """
    Generate a mask based on the darkest fraction of an image.

    .. versionadded:: 2.1.17

    Parameters
    ----------
    im : ndarray of shape (N,M)
        image used to generate a mask
    threshold: float, optional
        fraction of the lowest values to be masked, default = 10%

    Returns
    -------
    mask : boolean, ndarrays of shape (N,M)
        Mask that evaluates to True on valid pixels.
    """
    # Find the median of the highest intensity value of the image to avoid hot spots
    lower_limit = threshold * np.median(np.maximum(im, 0))
    # generate a mask
    return im > lower_limit
