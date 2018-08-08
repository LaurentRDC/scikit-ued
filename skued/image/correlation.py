# -*- coding: utf-8 -*-
"""
Image correlation and related functions
=======================================
"""
from functools import partial

import numpy as np
from scipy.fftpack import next_fast_len, fftn, ifftn
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

def mnxc2(arr1, arr2, m1 = None, m2 = None, mode = 'full', axes = (0, 1), out = None, overlap_ratio = 3/10):
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
    overlap_ratio : float, optional
        Maximum allowed overlap ratio between masks. The correlation at pixels with overlap ratio higher
        than this threshold will be zeroed.
        
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
    # This function is as close as possible to a direct translation of Dirk Padfield's
    # MATLAB implementation `normxcorr2_masked` function. You can find the source code here: 
    # http://www.dirkpadfield.com/Home/MaskedFFTRegistrationCode.zip
    # For backwards-compatibility reasons, zeroing of pixels based on mask overlap is done at the
    # end of this function, and not inside the translation_registration function.

    # TODO: implement for complex arrays
    # TODO: implement multidims

    if mode not in {'full', 'same'}:
        raise ValueError("Correlation mode {} is not valid.".format(mode))

    if len(axes) != 2:
        raise ValueError('`axes` parameter must be 2-tuple, not `{}`'.format(axes))

    fixed_image = np.array(arr1, dtype = np.float)
    fixed_mask = np.zeros_like(fixed_image, dtype = np.bool) if (m1 is None) else np.array(m1, dtype = np.bool)
    moving_image = np.array(arr2, dtype = np.float)
    moving_mask = np.zeros_like(moving_image, dtype = np.bool) if (m2 is None) else np.array(m2, dtype = np.bool)
    eps = np.finfo(np.float).eps

    # Note that by Padfield's implementation, we require that masks be 1 where pixels are
    # valid. Therefore, we must invert the bits.
    np.logical_not(fixed_mask, out = fixed_mask)
    np.logical_not(moving_mask, out = moving_mask)

    # Determine final size along transformation axes
    # Note that it might be faster to conmpute Fourier transform in a slightly larger shape (`fast_shape`)
    # Then, after all fourier transforms are done, we slice back to `final_shape` using `final_slice`.
    s1, s2 = tuple(fixed_image.shape[ax] for ax in axes), tuple(moving_image.shape[ax] for ax in axes)
    final_shape = tuple( ax1 + ax2 - 1 for ax1, ax2 in zip(s1, s2))
    fast_shape = tuple( map(next_fast_len, final_shape) )
    final_slice = tuple([slice(0, int(sz)) for sz in final_shape])

    fft = partial(rfft2, s = fast_shape, axes = axes, **FFTOPS)
    ifft = partial(irfft2, s = fast_shape, axes = axes, **FFTOPS)

    fixed_image[np.logical_not(fixed_mask)] = 0.0
    moving_image[np.logical_not(moving_mask)] = 0.0

    rotated_moving_image = np.rot90(moving_image, 2, axes = axes)
    rotated_moving_mask = np.rot90(moving_mask, 2, axes = axes)

    fixed_fft = fft(fixed_image)
    rotated_moving_fft = fft(rotated_moving_image)
    fixed_mask_fft = fft(fixed_mask)
    rotated_moving_mask_fft = fft(rotated_moving_mask)

    # Calculate overlap of masks at every point in the convolution
    # Locations with high overlap should not be taken into account.
    number_overlap_masked_px = ifft(rotated_moving_mask_fft * fixed_mask_fft)
    number_overlap_masked_px[:] = np.round(number_overlap_masked_px)
    number_overlap_masked_px[:] = np.fmax(number_overlap_masked_px, eps)
    masked_correlated_fixed_fft = ifft(rotated_moving_mask_fft * fixed_fft)
    masked_correlated_rotated_moving_fft = np.real(ifft(fixed_mask_fft * rotated_moving_fft))

    numerator = ifft(rotated_moving_fft * fixed_fft)
    numerator -= masked_correlated_fixed_fft * masked_correlated_rotated_moving_fft / number_overlap_masked_px

    fixed_squared_fft = fft(np.square(fixed_image))
    fixed_denom = ifft(rotated_moving_mask_fft * fixed_squared_fft)
    fixed_denom -= np.square(masked_correlated_fixed_fft) / number_overlap_masked_px
    fixed_denom[:] = np.fmax(fixed_denom, 0.0)

    rotated_moving_squared_fft = fft(np.square(rotated_moving_image))
    moving_denom = ifft(fixed_mask_fft * rotated_moving_squared_fft)
    moving_denom -= np.square(masked_correlated_rotated_moving_fft) / number_overlap_masked_px
    moving_denom[:] = np.fmax(moving_denom, 0.0)

    denom = np.sqrt(fixed_denom * moving_denom)

    # Slice back to expected convolution shape
    numerator = numerator[final_slice]
    denom = denom[final_slice]
    number_overlap_masked_px = number_overlap_masked_px[final_slice]

    if mode == 'same':
        denom = _centered(denom, fixed_image.shape, axes = axes)
        numerator = _centered(numerator, fixed_image.shape, axes = axes)
        number_overlap_masked_px = _centered(number_overlap_masked_px, fixed_image.shape, axes = axes)

    if out is None:
        out = np.zeros_like(denom)

    tol = 1e3 * eps * np.max(np.abs(denom), axis = axes, keepdims = True)
    nonzero_indices = denom > tol
    out[nonzero_indices] = numerator[nonzero_indices] / denom[nonzero_indices]
    np.clip(out, a_min = -1, a_max = 1, out = out)

    # Apply overlap ratio threshold
    number_px_threshold = overlap_ratio * np.max(number_overlap_masked_px, axis = axes, keepdims = True)
    out[number_overlap_masked_px < number_px_threshold] = 0.0 

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
