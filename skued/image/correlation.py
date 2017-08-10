"""
Image correlation and related functions
=======================================
"""
from functools import partial
import numpy as np
from scipy.fftpack import next_fast_len

try:
    from pyfftw.interfaces.numpy_fft import rfft2, irfft2
except ImportError:
    from numpy.fft import rfft2, irfft2

from ..array_utils import mirror

EPS = np.finfo(np.float).eps

def mnxc2(arr1, arr2, m1 = None, m2 = None, mode = 'full', axes = (0, 1)):
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
	# TODO: implement multidims

	if mode not in {'full', 'same'}:
		raise ValueError("Correlation mode {} is not valid.".format(mode))

	if len(axes) != 2:
		raise ValueError('`axes` parameter must be 2-tuple, not `{}`'.format(axes))

	arr1, arr2 = np.array(arr1, dtype = np.float), np.array(arr2, dtype = np.float)

	# Determine final size along transformation axes
	# TODO: compare with using next_fast_len and without
	s1, s2 = tuple(arr1.shape[ax] for ax in axes), tuple(arr2.shape[ax] for ax in axes)
	final_shape = tuple( next_fast_len(ax1 + ax2 - 1) for ax1, ax2 in zip(s1, s2))

	fft = partial(rfft2, s = final_shape, axes = axes)
	ifft = partial(irfft2, s = final_shape, axes = axes)
	rot180 = lambda arr : mirror(mirror(arr, axes[0]), axes[1]) 	# numpy.flip not available in numpy 1.11

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
	arr2[:] = rot180(arr2)
	m2[:] = rot180(m2)

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
	out[np.logical_or(out > 1, out < -1)] = 0

	if mode == 'full':
		return out
	elif mode == 'same':
		return _centered(out, arr1.shape, axes = axes)

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