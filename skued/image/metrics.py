# -*- coding: utf-8 -*-
"""
Image mask routines
===================
"""
from collections.abc import Iterable
from itertools import repeat

import numpy as np

# array_stream decorator ensures that input images are cast to ndarrays
from npstreams import array_stream, imean, istd, itercopy, last, peek, prod


@array_stream
def isnr(images, fill_value=0.0):
    """
    Streaming, pixelwise signal-to-noise ratio (SNR).

    Parameters
    ----------
    images : iterable of ndarray
        These images should represent identical measurements. ``images`` can also be a generator.
    fill_value : float, optional
        Division-by-zero results will be filled with this value.

    Yields
    ------
    snr : `~numpy.ndarray`
        Pixelwise signal-to-noise ratio

    See Also
    --------
    snr_from_collection : pixelwise signal-to-noise ratio from a collection of measurements
    """
    first, images = peek(images)
    snr = np.empty_like(first)

    images1, images2 = itercopy(images, 2)
    for mean, std in zip(imean(images1), istd(images2)):
        valid = std != 0
        snr[valid] = mean[valid] / std[valid]
        snr[np.logical_not(valid)] = fill_value
        yield snr


@array_stream
def snr_from_collection(images, fill_value=0.0):
    """
    Signal-to-noise ratio (SNR) on a per-pixel basis, for images in a collection.
    These images should represent identical measurements.

    SNR is defined as :math:`snr = \\mu/\sigma` where :math:`\\mu`
    is the average pixel value and :math:`\sigma` is the standard deviation of that pixel value.

    This function operates in constant-memory; it is therefore safe to use on a large collection
    of images (>10GB).

    Parameters
    ----------
    images : iterable of ndarray
        These images should represent identical measurements. ``images`` can also be a generator.
    fill_value : float, optional
        Division-by-zero results will be filled with this value.

    Returns
    -------
    snr : `~numpy.ndarray`
        Pixelwise signal-to-noise ratio

    See Also
    --------
    isnr : streaming signal-to-noise ratio
    """
    return last(isnr(images, fill_value=fill_value))


@array_stream
def mask_from_collection(images, px_thresh=(0, 3e4), std_thresh=None):
    """
    Determine binary mask from a set of images. These images should represent identical measurements, e.g. a set
    of diffraction patterns before photoexcitation. Pixels are rejected on the following two criteria:

        * Pixels with a value above a certain threshold or below zero, for any image in the set, are considered dead;
        * Pixels with a cumulative standard deviation above a certain threshold are considered uncertain.

    This function operates in constant-memory; it is therefore safe to use on a large collection
    of images (>10GB).

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
    reject all negative pixels only, set ``px_thresh = (0, numpy.inf)``.
    """
    if isinstance(px_thresh, Iterable):
        min_int, max_int = min(px_thresh), max(px_thresh)
    else:
        min_int, max_int = None, px_thresh

    first, images = peek(images)
    mask = np.zeros_like(first, dtype=bool)  # 0 = False

    if std_thresh is not None:
        images, images_for_std = itercopy(images)
        std_calc = istd(images_for_std)
    else:
        std_calc = repeat(np.inf)

    for image, std in zip(images, std_calc):

        mask[image > max_int] = True

        if std_thresh is not None:
            mask[std > std_thresh] = True

        if min_int is not None:
            mask[image < min_int] = True

    return mask


def combine_masks(*masks):
    """
    Combine multiple pixel masks into one. This assumes that pixel masks evaluate
    to ``True`` on valid pixels and ``False`` on invalid pixels.

    Returns
    -------
    combined : `~numpy.ndarray`, dtype bool
    """
    # By multiplying boolean arrays, values of False propagate
    return prod(masks, dtype=bool)


def mask_image(image, mask, fill_value=0, copy=True):
    """
    Fill invalid pixels in an image with another value, according to a pixel mask.
    While this function has simply functionality, it is ideal for integrating into a
    pipeline with ``functools.partial``.

    Parameters
    ---------
    image : `~numpy.ndarray`

    mask : `~numpy.ndarray`
        Boolean array. ``mask`` should evaluate to ``True`` on valid pixels.
    fill_value : float, optional
        Invalid pixels fill value.
    copy : bool, optional
        If True (default), ``image`` is copied before masking. If False, ``image`` is modified in-place.

    Returns
    -------
    masked : `~numpy.ndarray`
        Masked image. If ``copy = True``, masked points to the same object as ``image``.
    """
    image = np.array(image, copy=copy)
    mask = np.asarray(mask, dtype=bool)

    image[np.logical_not(mask)] = fill_value
    return image


def triml(array, percentile, axis=None, fill_value=0):
    """
    Trim values in an array that fall below (i.e. to the left) a certain percentile.

    Parameters
    ----------
    array : `~numpy.ndarray`
        Array to be trimmed, typically an image.
    percentile : float in range [0, 100]
        Percentile below which array elements are set to ``fill_value``.
    axis : int or None, optional
        Axis along which to trim data. If None (default), compute over the whole array.
    fill_value : float, optional
        Trimmed array elements are replaced with this value.

    Returns
    -------
    trimmed : `~numpy.ndarray`
        Trimmed array of the same shape as ``array``.

    See Also
    --------
    trimr : trim values in percentiles above a specific percentile.
    """
    array = np.array(array)
    val_percentile = np.percentile(array, q=float(percentile), axis=axis, keepdims=True)
    array[array < val_percentile] = fill_value
    return array


def trimr(array, percentile, axis=None, fill_value=0):
    """
    Trim values in an array that fall above (i.e. to the right) a certain percentile.

    Parameters
    ----------
    array : `~numpy.ndarray`
        Array to be trimmed, typically an image.
    percentile : float in range [0, 100]
        Percentile above which array elements are set to ``fill_value``.
    axis : int or None, optional
        Axis along which to trim data. If None (default), compute over the whole array.
    fill_value : float, optional
        Trimmed array elements are replaced with this value.

    Returns
    -------
    trimmed : `~numpy.ndarray`
        Trimmed array of the same shape as ``array``.

    See Also
    --------
    triml : trim values in percentiles below a specific percentile.
    """
    array = np.array(array)
    val_percentile = np.percentile(array, q=float(percentile), axis=axis, keepdims=True)
    array[array > val_percentile] = fill_value
    return array
