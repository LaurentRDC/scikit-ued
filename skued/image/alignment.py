# -*- coding: utf-8 -*-
"""
Module concerned with alignment of diffraction images
=====================================================
"""
from functools import partial

import numpy as np
from scipy.ndimage import shift as subpixel_shift
from skimage.feature import register_translation
from skimage.filters import gaussian

from npstreams import array_stream, peek

from .correlation import mnxc2

non = lambda s: s if s < 0 else None
mom = lambda s: max(0, s)
def shift_image(arr, shift, fill_value = 0):
    """ 
    Shift an image. Subpixel resolution shifts are also possible.

    Parameters
    ----------
    arr : `~numpy.ndarray`
        Array to be shifted.
    shift : array_like, shape (2,)
        Shifts in the x and y directions, respectively. Shifts can be of sub-pixel value,
        in which case interpolation is used.
    fill_value : numerical, optional
        Edges will be filled with `fill_value` after shifting. 

    Returns
    -------
    out : `~numpy.ndarray`
        Shifted array. The type of the shifted array will be the smallest size
        that accomodates the types of `arr` and `fill_value`.
    
    See Also
    --------
    scipy.ndimage.shift : shift an image via interpolation
    """
    # Since the fill value is often NaN, but arrays may be integers
    # We need to promote the final type to smallest coherent type
    final_type = np.promote_types(arr.dtype, np.dtype(type(fill_value)))
    output = np.full_like(arr, fill_value = fill_value, dtype = final_type)

    # Floating point shifts are much slower

    j, i = tuple(shift)
    if (int(i) != i) or (int(j) != j):	# shift is float
        subpixel_shift(arr, (i, j), output = output, cval = fill_value)
        return output
    
    i, j = int(i), int(j)

    dst_slices = [slice(None, None)] * arr.ndim
    src_slices = [slice(None, None)] * arr.ndim

    for s, ax in zip((i, j), (0, 1)):
        dst_slices[ax] = slice(mom(s), non(s))
        src_slices[ax] = slice(mom(-s), non(-s))

    output[dst_slices] = arr[src_slices]
    return output

@array_stream
def itrack_peak(images, row_slice = None, col_slice = None, precision = 1/10):
    """
    Generator function that tracks a diffraction peak in a stream of images.
    
    Parameters
    ----------
    images : iterable of array-like
        Iterable of diffraction images. This function also supports generators.
    row_slice : slice or None, optional
        Slice object for which image rows to use. If None (default), all rows are used. 
    col_slice : slice or None, optional
        Slice object for which image columns to use. If None (default), all columns are used.
    precision : float, optional
        Precision of the tracking in pixels. A precision of 1/10 (default) means that
        the tracking will be precise up to 1/10 of a pixel.
    
    Yields
    ------
    shift : `~numpy.ndarray`, shape (2,)
        [row, col] shifts of the peak with respect to the position of this peak
        in the first image.
    """
    if row_slice is None:
        row_slice = np.s_[:]
    
    if col_slice is None:
        col_slice = np.s_[:]

    first = next(images)
    
    # The shift between the first image and itself needs not
    # be computed!
    yield np.array((0.0, 0.0))

    ref = np.array(first[row_slice, col_slice], copy = True)
    sub = np.empty_like(ref)

    for image in images:
        sub[:] = image[row_slice, col_slice]
        shift, *_ = register_translation(ref, sub, upsample_factor = int(1/precision))
        yield np.asarray(shift)

def align(image, reference, mask = None, fill_value = 0.0, fast = True):
    """
    Align a diffraction image to a reference. Subpixel resolution available.

    Parameters
    ----------
    image : `~numpy.ndarray`, shape (M,N)
        Image to be aligned.
    reference : `~numpy.ndarray`, shape (M,N)
        `image` will be align onto the `reference` image.
    mask : `~numpy.ndarray` or None, optional
        Mask that evaluates to True on invalid pixels of the array `image`.
    fill_value : float, optional
        Edges will be filled with `fill_value` after alignment.
    fast : bool, optional
        If True (default), alignment is done on images cropped to half
        (one quarter area). Disable for small images, e.g. 256x256.
    
    Returns
    -------
    aligned : `~numpy.ndarray`, shape (M,N)
        Aligned image.
    
    See Also
    --------
    ialign : generator of aligned images
    """
    shift = diff_register(image, reference = reference, mask = mask, crop = fast)
    return shift_image(image, shift, fill_value = fill_value)

@array_stream
def ialign(images, reference = None, mask = None, fill_value = 0.0, fast = True):
    """
    Generator of aligned diffraction images.

    Parameters
    ----------
    images : iterable
        Iterable of ndarrays of shape (N,M)
    reference : `~numpy.ndarray`, shape (M,N)
        Images in `images` will be aligned onto the `reference` image. If
        'reference' is None (default), the first image in the 'images' stream
        is used as a reference
    mask : `~numpy.ndarray` or None, optional
        Mask that evaluates to True on invalid pixels.
    fill_value : float, optional
        Edges will be filled with `fill_value` after alignment.
    fast : bool, optional
        If True (default), alignment is done on images cropped to half
        (one quarter area). Disable for small images, e.g. 256x256.

    Yields
    ------
    aligned : `~numpy.ndarray`
        Aligned image. If `reference` is None, the first aligned image is the reference.

    See Also
    --------
    skued.align : align a single diffraction pattern onto a reference.
    """
    images = iter(images)
    
    if reference is None:
        reference = next(images)
        yield reference

    yield from map(partial(align, reference = reference, mask = mask, fill_value =  fill_value, fast = fast), images)


def _crop_to_half(image, copy = False):
    nrows, ncols = np.array(image.shape)/4
    return np.array(image[int(nrows):-int(nrows), int(ncols):-int(ncols)], copy = copy)

# TODO: add option to upsample, akin to skimage.feature.register_translation
#		Could this be done initially by zero-padding?
#		See https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/register_translation.py#L109
def diff_register(image, reference, mask = None, crop = True, sigma = 5):
    """
    Register translation of diffraction patterns by masked 
    normalized cross-correlation.
    
    Parameters
    ----------
    image : iterable
        Iterable of ndarrays of shape (N,M)
    reference : `~numpy.ndarray`
        This is the reference image to which `image` will be aligned. 
    mask : `~numpy.ndarray` or None, optional
        Mask that evaluates to True on invalid pixels of the array `image`.
    crop : bool, optional
        If True (default), ``image`` and ``reference`` are cropped to one
        quarter of their areas; this results in faster execution at the expense of 
        precision. Disable for small images.
    sigma : float or None, optional
        Standard deviation for Gaussian kernel with which to smooth 
        ``image`` and ``reference``. If None, no smoothing is performed.
    
    Returns
    -------
    shift : `~numpy.ndarray`, shape (2,), dtype float
        Shift in rows and columns. The ordering is compatible with :func:`shift_image`
    
    References
    ----------
    .. [PADF] Dirk Padfield. Masked Object Registration in the Fourier Domain. 
        IEEE Transactions on Image Processing, vol.21(5), pp. 2706-2718, 2012. 
    """
    if mask is None:
        mask = np.zeros_like(image, dtype = np.bool)
    
    if crop:
        image = _crop_to_half(image, copy = True)
        reference = _crop_to_half(reference, copy = True)
        mask = _crop_to_half(mask, copy = True)

    # Diffraction images register better with some filtering
    if sigma:
        image = gaussian(image, sigma, preserve_range = True)
        reference = gaussian(reference, sigma, preserve_range = True)

    # Contrary to Padfield, we do not have to crop out the edge
    # since we are using the 'valid' correlation mode.
    xcorr = mnxc2(reference, image, mask, mode = 'same')

    # Generalize to the average of multiple maxima
    maxima = np.transpose(np.nonzero(xcorr == xcorr.max()))
    center = np.mean(maxima, axis = 0)
    
    # Due to centering of mnxc2, +1 is required
    # TODO: was this due to wrong output shape
    # 		of mnxc2?
    shift_row_col = center - np.array(xcorr.shape)/2 + 1
    return -shift_row_col[::-1]	# Reversing to be compatible with shift_image