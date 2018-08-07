# -*- coding: utf-8 -*-
"""
Image correlation and related functions
=======================================
"""
from functools import partial

import numpy as np
from scipy.fftpack import next_fast_len, fftn, ifftn, fft2, ifft2
from scipy.signal import fftconvolve

from ..array_utils import mirror

FFTOPS = {}
try:
    from pyfftw.interfaces.numpy_fft import rfft2, irfft2
    FFTOPS['threads'] = 2
except ImportError:
    from numpy.fft import rfft2, irfft2


EPS = np.finfo(np.float).eps

def xcorr(arr1, arr2, mode = 'full', axes = None):
    """ 
    Cross-correlation between two N-dimensional arrays.
    Support for cross-correlation along specific axes as well.

    Parameters
    ----------
    arr1 : `~numpy.ndarray`
        First array. 
    arr2 : `~numpy.ndarray`
        Second array. It should have the same number of dimensions
        as ``arr1``.
    mode : {'full', 'same'}, optional
        'full':
            By default, mode is 'full'.  This returns the convolution
            at each point of overlap, with an output shape of (N+M-1,M+N-1). At
            the end-points of the convolution, the signals do not overlap
            completely, and boundary effects may be seen.
        'same':
            Mode 'same' returns output of length ``max(M, N)``. Boundary
            effects are still visible.
    axes : None or tuple of ints, optional
        Axes over which to compute the cross-correlation. If None,
        cross-correlation is computed over all axes.

    Returns
    -------
    out : `~numpy.ndarray`
        Cross-correlation. If ``arr1`` or ``arr2`` is complex, then
        ``out`` will be complex as well; otherwise, ``out`` will be
        of float type.

    Notes
    -----
    If ``pyfftw`` is installed, its routines will be preferred.
    
    Raises
    ------
    ValueError : mode argument is invalid
    ValueError : if either ``arr1`` or ``arr2`` is complex.

    See Also
    --------
    mnxc2 : masked normalized cross-correlation between images.
    """
    if mode not in {'full', 'same'}:
        raise ValueError('Unexpected cross-correlation mode {}'.format(mode))
    
    if axes is None:
        axes = tuple(range(arr1.ndim))

    arr1, arr2 = np.asarray(arr1), np.asarray(arr2)

    # Determine final size along transformation axes
    # To speed up FFT, shape of Fourier transform might be slightly larger
    # then slice back before returning
    s1 = tuple(arr1.shape[ax] for ax in axes)
    s2 = tuple(arr2.shape[ax] for ax in axes)
    final_shape = tuple( ax1 + ax2 - 1 for ax1, ax2 in zip(s1, s2))
    fast_shape = tuple(map(next_fast_len, final_shape))
    final_slice = tuple([slice(0, int(sz)) for sz in final_shape])

    F1 = fftn(arr1, shape = fast_shape, axes = axes)
    F2 = fftn(np.conj(mirror(arr2, axes = axes)), shape = fast_shape, axes = axes)
    xc = ifftn(F1 * F2)[final_slice]

    if mode == 'same':
        return _centered(xc, arr1.shape, axes = axes)
    else:
        return xc

def mnxc2(arr1, arr2, m1 = None, m2 = None, mode = 'full', axes = (0, 1), out = None):
    """
    Masked normalized cross-correlation (MNXC) between two images or stacks of images.

    Parameters
    ----------
    arr1 : `~numpy.ndarray`, shape (M,N)
        Reference, or 'fixed-image' in the language of _[PADF]. This array can also
        be a stack of images; in this case, the cross-correlation
        is computed along the two axes passed to ``axes``.
    arr2 : `~numpy.ndarray`, shape (M,N)
        Moving image. This array can also be a stack of images; 
        in this case, the cross-correlation is computed along the 
        two axes passed to ``axes``.
    m1 : `~numpy.ndarray`, shape (M,N) or None, optional
        Mask of `arr1`. The mask should evaluate to `True`
        (or 1) on invalid pixels. If None (default), no mask
        is used.
    m2 : `~numpy.ndarray`, shape (M,N) or None, optional
        Mask of `arr2`. The mask should evaluate to `True`
        (or 1) on invalid pixels. If None (default), `m2` is 
        taken to be the same as `m1`.	
    mode : {'full', 'same'}, optional
        'full':
            By default, mode is 'full'.  This returns the convolution
            at each point of overlap, with an output shape of (N+M-1,M+N-1). At
            the end-points of the convolution, the signals do not overlap
            completely, and boundary effects may be seen.
        'same':
            Mode 'same' returns output of length ``max(M, N)``. Boundary
            effects are still visible.
    axes : 2-tuple of ints, optional
        Axes along which to compute the cross-correlation.
    out : `~numpy.ndarray` or None, optional
        If not None, the results will be stored in `out`. If None, a new array
        is returned.
        
    Returns
    -------
    out : `~numpy.ndarray`
        Masked, normalized cross-correlation. If images are real-valued, then `out` will be
        real-valued as well. For complex input, `out` will be complex as well.
    
    See Also
    --------
    xcorr : cross-correlation between two N-dimensional arrays.
        
    References
    ----------
    .. [PADF] Dirk Padfield. Masked Object Registration in the Fourier Domain. 
        IEEE Transactions on Image Processing, vol.21(5), pp. 2706-2718 (2012). 
    """
    # TODO: implement for complex arrays
    # TODO: implement multidims

    if mode not in {'full', 'same'}:
        raise ValueError("Correlation mode {} is not valid.".format(mode))

    if len(axes) != 2:
        raise ValueError('`axes` parameter must be 2-tuple, not `{}`'.format(axes))

    arr1, arr2 = np.array(arr1, dtype = np.float), np.array(arr2, dtype = np.float)

    # Determine final size along transformation axes
    # TODO: compare with using next_fast_len and without
    s1, s2 = tuple(arr1.shape[ax] for ax in axes), tuple(arr2.shape[ax] for ax in axes)
    final_shape = tuple( ax1 + ax2 - 1 for ax1, ax2 in zip(s1, s2))
    fast_shape = tuple( map(next_fast_len, final_shape) )
    final_slice = tuple([slice(0, int(sz)) for sz in final_shape])

    fft = partial(rfft2, s = fast_shape, axes = axes, **FFTOPS)
    ifft = partial(irfft2, s = fast_shape, axes = axes, **FFTOPS)

    if m1 is None:
        m1 = np.zeros_like(arr1, dtype = np.bool)
    else:
        m1 = np.array(m1, dtype = np.bool)

    if m2 is None:
        m2 = np.array(m1, dtype = np.bool)
    else:
        m2 = np.array(m2, dtype = np.bool)

    arr1[m1] = 0.0
    arr2[m2] = 0.0

    # Rotation in real-space instead of conjugation in fourier domain
    # because we might be using rfft instead of complex fft
    arr2[:] = np.rot90(arr2, k = 2)
    m2[:] = np.rot90(m2, k = 2)

    F1 = fft(arr1)
    F2s = fft(arr2)

    M1 = fft(np.logical_not(m1))
    M2s = fft(np.logical_not(m2))

    iM1M2s = ifft(M1 * M2s)
    iM1M2s[:] = np.rint(iM1M2s)
    iM1M2s[:] = np.maximum(iM1M2s, EPS)

    iF1M2s = ifft(F1 * M2s)
    iM1F2s = ifft(M1 * F2s)

    # I have noticed no clear performance boost by storing
    # repeated calculation (e.g. ifft(M1 * M2s)); however, the following
    # is already hard enough to read...
    numerator = ifft(F1 * F2s)
    numerator -= iF1M2s * iM1F2s / iM1M2s

    denominator = ifft(fft(np.square(arr1)) * M2s) - iF1M2s**2/iM1M2s
    denominator *= ifft(M1*fft(np.square(arr2))) - iM1F2s**2/iM1M2s
    denominator[:] = np.clip(denominator, a_min = 0, a_max = None)
    denominator[:] = np.sqrt(denominator)

    # Slice back to convolution shape
    numerator = numerator[final_slice]
    denominator = denominator[final_slice]

    if mode == 'same':
        denominator = _centered(denominator, arr1.shape, axes = axes)
        numerator = _centered(numerator, arr1.shape, axes = axes)

    if out is None:
        out = np.zeros_like(denominator)

    nonzero = np.nonzero(denominator)
    out[nonzero] = numerator[nonzero] / denominator[nonzero]
    out[np.logical_or(out > 1, out < -1)] = 0

    return out

def _centered(arr, newshape, axes = (0, 1)):
	# Return the center newshape portion of the array.
	newshape = np.asarray(newshape)
	currshape = np.array(arr.shape)

	slices = [slice(None, None)] * arr.ndim

	for ax in axes:
		startind = (currshape[ax] - newshape[ax]) // 2
		endind = startind + newshape[ax]
		slices[ax] = slice(startind, endind)

	return arr[tuple(slices)]

def normxcorr2_masked(fixed_image, moving_image, fixed_mask, moving_mask):
    """
    Masked normalized cross-correlation (MNXC) between two images or stacks of images.

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
    moving_mask : `~numpy.ndarray`, shape (M,N)
        Mask of `moving_image`. The mask should evaluate to `True`
        (or 1) on valid pixels. 
        
    Returns
    -------
    out : `~numpy.ndarray`
        Masked, normalized cross-correlation. If images are real-valued, then `out` will be
        real-valued. For complex input, `out` will be complex as well.
    masked_overlap : `~numpy.ndarray`
        
    References
    ----------
    .. [PADF] Dirk Padfield. Masked Object Registration in the Fourier Domain. 
        IEEE Transactions on Image Processing, vol.21(5), pp. 2706-2718 (2012). 
    """
    # TODO: support axes
    #       problem is only the final resizing
    axes = (0, 1)
    eps = np.finfo(np.float).eps

    fixed_image = np.array(fixed_image, dtype = np.float)
    fixed_mask = np.array(fixed_mask, dtype = np.bool)
    moving_image = np.array(moving_image, dtype = np.float)
    moving_mask = np.array(moving_mask, dtype = np.bool)

    # Determine final size along transformation axes
    s1, s2 = tuple(fixed_image.shape[ax] for ax in axes), tuple(moving_image.shape[ax] for ax in axes)
    final_shape = tuple( ax1 + ax2 - 1 for ax1, ax2 in zip(s1, s2))
    fast_shape = tuple( map(next_fast_len, final_shape) )
    final_slice = tuple([slice(0, int(sz)) for sz in final_shape])

    fft = partial(fft2, shape = fast_shape, axes = axes)
    ifft = partial(ifft2, shape = fast_shape, axes = axes)

    fixed_image[np.logical_not(fixed_mask)] = 0.0
    moving_image[np.logical_not(moving_mask)] = 0.0

    rotated_moving_image = np.rot90(moving_image, 2)
    rotated_moving_mask = np.rot90(moving_mask, 2)

    fixed_fft = fft(fixed_image)
    rotated_moving_fft = fft(rotated_moving_image)
    fixed_mask_fft = fft(fixed_mask)
    rotated_moving_mask_fft = fft(rotated_moving_mask)

    number_overlap_masked_px = np.real(ifft(rotated_moving_mask_fft * fixed_mask_fft))
    number_overlap_masked_px[:] = np.round(number_overlap_masked_px)
    number_overlap_masked_px[:] = np.maximum(number_overlap_masked_px, eps)
    masked_correlated_fixed_fft = np.real(ifft(rotated_moving_mask_fft * fixed_fft))
    masked_correlated_rotated_moving_fft = np.real(ifft(fixed_mask_fft * rotated_moving_fft))

    numerator = np.real(ifft(rotated_moving_fft * fixed_fft)) 
    numerator -= masked_correlated_fixed_fft * masked_correlated_rotated_moving_fft / number_overlap_masked_px

    fixed_squared_fft = fft(np.square(fixed_image))
    fixed_denom = np.real(ifft(rotated_moving_mask_fft * fixed_squared_fft))
    fixed_denom -= np.square(masked_correlated_fixed_fft) / number_overlap_masked_px
    fixed_denom[:] = np.maximum(fixed_denom, 0.0)

    rotated_moving_squared_fft = fft(np.square(rotated_moving_image))
    moving_denom = np.real(ifft(fixed_mask_fft * rotated_moving_squared_fft))
    moving_denom -= np.square(masked_correlated_rotated_moving_fft) / number_overlap_masked_px
    moving_denom[:] = np.maximum(moving_denom, 0.0)

    denom = np.sqrt(fixed_denom * moving_denom)

    out = np.zeros_like(numerator)
    tol = 1e3 * eps * np.max(np.abs(denom))
    nonzero_indices = denom > tol
    out[nonzero_indices] = numerator[nonzero_indices] / denom[nonzero_indices]
    np.clip(out, a_min = -1, a_max = 1, out = out)
    
    return (np.resize(out, final_shape), 
            np.resize(number_overlap_masked_px, final_shape))