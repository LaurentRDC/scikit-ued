# -*- coding: utf-8 -*-
"""
Module concerned with alignment of diffraction images
=====================================================
"""
from itertools import product
import numpy as np
from skimage.feature import register_translation
from skimage.filters import gaussian
from warnings import warn
from .correlation import mnxc2

non = lambda s: s if s < 0 else None
mom = lambda s: max(0, s)
def shift_image(arr, shift, fill_value = 0):
	""" 
	Shift an image at a 1-pixel resolution. Shift

	Parameters
	----------
	arr : `~numpy.ndarray`
		Array to be shifted.
	shift : array_like, shape (2,)
		Shifts in the 
	fill_value : numerical, optional
		Edges will be filled with `fill_value` after shifting.

	Returns
	-------
	out : ndarray
	"""
	x, y = tuple(shift)
	x, y = int(round(x)), int(round(y))

	shifted = np.full_like(arr, fill_value = fill_value)
	shifted[mom(y):non(y), mom(x):non(x)] = arr[mom(-y):non(-y), mom(-x):non(-x)]
	return shifted

def align(image, reference, fill_value = 0.0):
	"""
	Align a diffraction image to a reference.

	Parameters
	----------
	image : `~numpy.ndarray`, shape (M,N)
		Image to be aligned.
	reference : `~numpy.ndarray`, shape (M,N)
		`image` will be align onto the `reference` image.
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
	shifts, *_ = register_translation(target_image = _crop_to_half(image), 
									  src_image = _crop_to_half(reference), 
									  upsample_factor = 4, space = 'real')
	return shift_image(image, -shifts, fill_value = fill_value)

def _crop_to_half(image, copy = False):
	nrows, ncols = np.array(image.shape)/4
	return np.array(image[int(nrows):-int(nrows), int(ncols):-int(ncols)], copy = copy)

def diff_register(image, reference, mask = None):
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
	
	Returns
	-------
	shift : `~numpy.ndarray`, shape (2,)
		Shift in rows and columns. Compatible with `shift_image`
	
	References
	----------
	.. [PADF] Dirk Padfield. Masked Object Registration in the Fourier Domain. 
		IEEE Transactions on Image Processing, vol.21(5), pp. 2706-2718, 2012. 
	"""
	if mask is None:
		mask = np.zeros_like(image, dtype = np.bool)

	cropped = _crop_to_half(image, copy = True)
	cropped_ref = _crop_to_half(reference, copy = True)
	cropped_mask = _crop_to_half(mask, copy = True)

	# Diffraction images register better with some filtering
	cropped[:] = gaussian(cropped, 5, preserve_range = True)
	cropped_ref[:] = gaussian(cropped_ref, 5, preserve_range = True)

	# Contrary to Padfield, we do not have to crop out the edge
	# since we are using the 'valid' correlation mode.
	xcorr = mnxc2(cropped_ref, cropped, cropped_mask, mode = 'same')

	# Generalize to the average of multiple maxima
	maxima = np.transpose(np.nonzero(xcorr == xcorr.max()))
	center = np.mean(maxima, axis = 0)
	
	# Due to centering of mnxc2, +1 is required
	shift_row_col = center - np.array(xcorr.shape)/2  + 1
	return -shift_row_col[::-1].astype(np.int)	# Reversing to be compatible with shift_image