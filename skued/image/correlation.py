"""
Image correlation and related functions
=======================================
"""
from functools import partial
import numpy as np
from numpy.fft import rfft2, irfft2, fft2, ifft2, ifftshift, fftshift
from scipy.fftpack import next_fast_len

EPS = max(np.finfo(np.float).eps, np.finfo(np.complex).eps)

def _crop_to_half(image):
	nrows, ncols = np.array(image.shape)/4
	return image[int(nrows):-int(nrows), int(ncols):-int(ncols)]

def mnxc2(arr1, arr2, m1 = None, m2 = None):
	"""
	Masked normalized cross-correlation (MNXC) between two images.

	Parameters
	----------
	arr1 : `~numpy.ndarray`, shape (M,N)
		Reference, or 'fixed-image' in the language of _[PADF].
	arr2 : `~numpy.ndarray`, shape (M,N)
		Moving image
	m1 : `~numpy.ndarray`, shape (M,N) or None, optional
		Mask of `arr1`. The mask should evaluate to `True`
		(or 1) on invalid pixels. If None (default), no mask
		is used.
	m2 : `~numpy.ndarray`, shape (M,N) or None, optional
		Mask of `arr2`. The mask should evaluate to `True`
		(or 1) on invalid pixels. If None (default), `m2` is 
		taken to be the same as `m1`.
		
	Returns
	-------
	out : `~numpy.ndarray`
		Masked, normalized cross-correlation. If images are real-valued, then `out` will be
		real-valued as well. For complex input, `out` will be complex as well.
		
	References
	----------
	.. [PADF] Dirk Padfield. Masked Object Registration in the Fourier Domain. 
		IEEE Transactions on Image Processing, vol.21(5), pp. 2706-2718, 2012. 
	"""
	# TODO: implement for complex arrays
	# TODO: implement over axes
	# TODO: implement multidims

	arr1, arr2 = np.array(arr1), np.array(arr2)

	# Determine final size along transformation axes
	# TODO: compare with using next_fast_len and without
	final_shape = tuple( next_fast_len(ax1 + ax2 - 1) for ax1, ax2 in zip(arr1.shape, arr2.shape))
	fft = partial(rfft2, s = final_shape)
	ifft = partial(irfft2, s = final_shape)
	
	if m1 is None:
		m1 = np.zeros_like(arr1, dtype = np.bool)
	else:
		m1 = np.array(m1)
	
	if m2 is None:
		m2 = np.array(m1)
	else:
		m2 = np.array(m2)
	
	arr1[m1] = 0.0
	arr2[m2] = 0.0

	# Rotation in real-space instead of conjugation in fourier domain
	# because we might be using rfft instead of complex fft
	arr2[:] = np.rot90(arr2, k = 2)
	m2[:] = np.rot90(m2, k = 2)
	
	F1 = fft(arr1)
	F2s = fft(arr2)

	M1 = fft(~m1)
	M2s = fft(~m2)

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

	denominator = ifft(fft(arr1*arr1) * M2s) - iF1M2s**2/iM1M2s
	denominator *= ifft(M1*fft(arr2*arr2)) - iM1F2s**2/iM1M2s
	denominator[:] = np.clip(denominator, a_min = 0, a_max = None)
	denominator[:] = np.sqrt(denominator)

	out = np.zeros_like(denominator)
	nonzero = np.nonzero(denominator)
	out[nonzero] = numerator[nonzero] / denominator[nonzero]

	out = _centered(out, arr1.shape)
	out[np.logical_or(out > 1, out < -1)] = 0
	return out

def register_translation(fixed_image, moving_image, fixed_mask = None, moving_mask = None):
	"""
	Determine the translation between two images using the masked normalized cross-correlation.

	Parameters
	----------
	fixed_image : `~numpy.ndarray`, shape (M,N)
		Reference image
	moving_image : `~numpy.ndarray`, shape (M,N)
		Moving image
	fixed_mask : `~numpy.ndarray`, shape (M,N) or None, optional
		Mask of `fixed_image`. The mask should evaluate to `True`
		(or 1) on invalid pixels. If None (default), no mask
		is used.
	moving_mask : `~numpy.ndarray`, shape (M,N) or None, optional
		Mask of `moving_image`. The mask should evaluate to `True`
		(or 1) on invalid pixels. If None (default), `moving_mask` is 
		taken to be the same as `fixed_mask`.
	
	Returns
	-------
	shift : `~numpy.ndarray`, shape (2,)
		Shift in [row, column]
	"""
	xcorr = mnxc2(fixed_image, moving_image, fixed_mask, moving_mask)

	# Generalize to the average of multiple maxima
	maxima = np.transpose(np.nonzero(xcorr == xcorr.max()))
	center = np.mean(maxima, axis = 0)
	
	# Due to centering of mnxc2, -1 is required
	shift_row_col = center - np.array(xcorr.shape)/2  + 1
	return shift_row_col[::-1]	# Reversing to be compatible with shift_image

def _centered(arr, newshape):
    # Return the center newshape portion of the array.
    newshape = np.asarray(newshape)
    currshape = np.array(arr.shape)
    startind = (currshape - newshape) // 2
    endind = startind + newshape
    myslice = [slice(startind[k], endind[k]) for k in range(len(endind))]
    return arr[tuple(myslice)]