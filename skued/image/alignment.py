# -*- coding: utf-8 -*-
"""
Module concerned with alignment of diffraction images
=====================================================
"""
from functools import partial

from os import cpu_count
import numpy as np
from warnings import warn
from npstreams import array_stream
from scipy import ndimage as ndi
from ..fft import with_skued_fft
from skimage.registration import phase_cross_correlation


@array_stream
@with_skued_fft
def itrack_peak(images, row_slice=None, col_slice=None, precision=1 / 10):
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

    ref = np.array(first[row_slice, col_slice], copy=True)
    sub = np.empty_like(ref)

    # scikit-image will use scipy.fft module
    # so we can increase the number of FFT workers
    # to get a performance speedup (50% in my tests)
    for image in images:
        sub[:] = image[row_slice, col_slice]

        shift = with_skued_fft(phase_cross_correlation)(
            reference_image=ref,
            moving_image=sub,
            return_error=False,
            upsample_factor=int(1 / precision),
        )
        yield np.asarray(shift)


def align(image, reference, mask=None, fill_value=0.0):
    """
    Align a diffraction image to a reference. Subpixel resolution available.

    Parameters
    ----------
    image : `~numpy.ndarray`, shape (M,N)
        Image to be aligned.
    reference : `~numpy.ndarray`, shape (M,N)
        `image` will be align onto the `reference` image.
    mask : `~numpy.ndarray` or None, optional
        Mask that evaluates to True on valid pixels of the array `image`.
    fill_value : float, optional
        Edges will be filled with `fill_value` after alignment.

    Returns
    -------
    aligned : `~numpy.ndarray`, shape (M,N)
        Aligned image.

    See Also
    --------
    ialign : generator of aligned images
    """
    shift = with_skued_fft(phase_cross_correlation)(
        reference_image=reference,
        moving_image=image,
        reference_mask=mask,
        return_error=False,
    )
    return ndi.shift(image, shift=shift, order=2, mode="constant", cval=fill_value)


@array_stream
def ialign(images, reference=None, mask=None, fill_value=0.0):
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
        Mask that evaluates to True on valid pixels.
    fill_value : float, optional
        Edges will be filled with `fill_value` after alignment.

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

    yield from map(
        partial(align, reference=reference, mask=mask, fill_value=fill_value),
        images,
    )
