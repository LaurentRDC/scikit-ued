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
from warnings import warn

from npstreams import array_stream, peek

from .correlation import mnxc2, mxcorr

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
        # Image must not be float16
        # because subpixel shifting involves interpolation
        subpixel_shift(arr.astype(np.float), (i, j), output = output, cval = fill_value)
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

def align(image, reference, mask = None, fill_value = 0.0, **kwargs):
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
    kwargs
        Keyword-arguments are passed to `skued.diff_register`.
    
    Returns
    -------
    aligned : `~numpy.ndarray`, shape (M,N)
        Aligned image.
    
    See Also
    --------
    ialign : generator of aligned images
    """
    if ('fast' in kwargs):
        warn('`fast` keyword arguments in `align` are deprecated. They will be ignored.', DeprecationWarning)
        kwargs.pop('fast', None)

    shift = diff_register(image, reference = reference, mask = mask, **kwargs)
    return shift_image(image, shift, fill_value = fill_value)

@array_stream
def ialign(images, reference = None, *args, **kwargs):
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
    kwargs
        Keyword-arguments are passed to `skued.diff_register`.

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

    yield from map(partial(align, reference = reference, *args, **kwargs), images)

# TODO: add option to upsample, akin to skimage.feature.register_translation
#		Could this be done initially by zero-padding?
#		See https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/register_translation.py#L109
def diff_register(image, reference, mask = None, **kwargs):
    """
    Register translation of diffraction patterns. Registration is either done using normalized cross-correlation
    or by masked normalized cross-correlation.
    
    Parameters
    ----------
    image : iterable
        Iterable of ndarrays of shape (N,M)
    reference : `~numpy.ndarray`
        This is the reference image to which `image` will be aligned. 
    mask : `~numpy.ndarray` or None, optional
        Mask that evaluates to True on invalid pixels of the array `image`.
    kwargs
        Keyword arguments are passed either to `skimage.feature.register_translation`, or 
        `skued.masked_register_translation` depending if a mask has been provided.
    
    Returns
    -------
    shift : `~numpy.ndarray`, shape (2,), dtype float
        Shift in rows and columns. The ordering is compatible with :func:`shift_image`
    
    References
    ----------
    .. [PADF] Dirk Padfield. Masked Object Registration in the Fourier Domain. 
        IEEE Transactions on Image Processing, vol.21(5), pp. 2706-2718, 2012. 
    """
    # Deprecation warning for unused keyword arguments
    if ('crop' in kwargs) or ('sigma' in kwargs):
        warn('`crop` and `sigma` keyword arguments in `diff_register` are deprecated. They will be ignored.', DeprecationWarning)
    kwargs.pop('crop', None); kwargs.pop('sigma', None)

    if mask is None:
        shifts, *_ = register_translation(image, reference, **kwargs)
    else:
        shifts = masked_register_translation(image, reference, 
                                             fixed_mask = np.logical_not(mask), 
                                             moving_mask = None, 
                                             **kwargs)
        
    return shifts

def masked_register_translation(fixed_image, moving_image, fixed_mask, moving_mask = None, overlap_ratio = 3/10):
    """
    Masked image translation registration.

    Parameters
    ----------
    fixed_image : `~numpy.ndarray`, shape (M,N)
        Reference, or 'fixed-image' in the language of _[PADF]. This array can also
        be a stack of images; in this case, the cross-correlation
        is computed along the two axes passed to ``axes``.
    moving_image : `~numpy.ndarray`, shape (M,N)
        Moving image. This array can also be a stack of images; 
        in this case, the cross-correlation is computed along the 
        two axes passed to ``axes``.
    fixed_mask : `~numpy.ndarray`, shape (M,N)
        Mask of `fixed_image`. The mask should evaluate to `True`
        (or 1) on valid pixels. 
    moving_mask : `~numpy.ndarray`, shape (M,N) or None, optional
        Mask of `moving_image`. The mask should evaluate to `True`
        (or 1) on valid pixels. If `None`, `fixed_mask` will be used in 
        place of `moving_mask`.
    overlap_ratio : float, optional
        TODO
        
    Returns
    -------
    shift : `~numpy.ndarray`, shape (2,), dtype float
        Shift in rows and columns. The ordering is compatible with :func:`shift_image`
    
    See Also
    --------
    skimage.feature.register_translation : image translation registration without masks.
        
    References
    ----------
    .. [PADF] Dirk Padfield. Masked Object Registration in the Fourier Domain. 
        IEEE Transactions on Image Processing, vol.21(5), pp. 2706-2718 (2012). 
    """
    fixed_mask = np.array(fixed_mask, dtype = np.bool)

    if moving_mask is None:
        moving_mask = np.array(fixed_mask)

    corr, overlap = mxcorr(fixed_image, moving_image, fixed_mask, moving_mask)

    number_px_threshold = overlap_ratio * np.max(overlap)
    corr[overlap < number_px_threshold] = 0.0

    # Generalize to the average of multiple maxima
    maxima = np.transpose(np.nonzero(corr == corr.max()))
    center = np.mean(maxima, axis = 0)
    shift = center - np.array(moving_image.shape) + 1
    return -shift