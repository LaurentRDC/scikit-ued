# -*- coding: utf-8 -*-
"""
Module concerned with alignment of diffraction images
"""
import numpy as np
from skimage.feature import register_translation as skimage_register_translation
import scipy.fftpack as fftpack

non = lambda s: s if s < 0 else None
mom = lambda s: max(0, s)
def shift_image(arr, shift, fill_value = 0):
	""" 
	Shift an image at a 1-pixel resolution.

	Parameters
	----------
	arr : ndarray
		Array to be shifted.
	shift : array_like, shape (2,)
		Shifts in x and y, i.e. (column shift, row shift)
	fill_value : numerical, optional
		Edges will be filled with `fill_value` after shifting.

	Returns
	-------
	out : ndarray
	"""
	y, x = tuple(shift)	# flip row/column
	x, y = int(round(x)), int(round(y))

	shifted = np.full_like(arr, fill_value = fill_value)
	shifted[mom(y):non(y), mom(x):non(x)] = arr[mom(-y):non(-y), mom(-x):non(-x)]
	return shifted

def align(images, reference = None, fill_value = 0.0, mask = None):
	"""
	Generator of aligned diffraction images.

	Parameters
	----------
	images : iterable
		Iterable of ndarrays of shape (N,M)
	reference : ndarray or None, optional
		If not None, this is the reference image to which all images will be aligned. Otherwise,
		images will be aligned to the first element of the iterable 'images'. 
	fill_value : numerical, optional
		Edges will be filled with `fill_value` after shifting.
	mask : `~numpy.ndarray` or None, optional
		The masks should evaluate to `False` (or 0) on invalid pixels.
    
	Yields
	------
	aligned : ndarray, ndim 2
		Aligned image

	Notes
	-----
	Diffraction images exhibit high symmetry in most cases, therefore images
	are cropped to a quarter of their size before alignment.
	"""
	# TODO: check if `images` is an array; if so, return the alignment (not yield)
	images = iter(images)
	
	if reference is None:
		reference = next(images)
		yield reference

	cropped_ref = _crop_to_half(reference)
	if mask is not None:
		cropped_mask = _crop_to_half(mask)
	else:
		cropped_mask = None

	def _align_im(image):
		# TODO: raise warning if error is too large?
		shifts = register_translation(target_image = _crop_to_half(image), 
								  	  src_image = cropped_ref,
									  src_mask = cropped_mask)
		return shift_image(image, -shifts, fill_value = fill_value)
		
	yield from map(_align_im, images)

def _crop_to_half(image):
	nrows, ncols = np.array(image.shape)/4
	return image[int(nrows):-int(nrows), int(ncols):-int(ncols)]

def masked_reg_trans(target_image, src_image, mask = None):
	""" 
	Image translation registration by masked cross-correlation, inspired from 
	`scikit-image`.

    Parameters
    ----------
    src_image : `~numpy.ndarray`
        Reference image.
    target_image : `~numpy.ndarray`
        Image to register.
	mask : `~numpy.ndarray` or None, optional
		The masks should evaluate to `False` (or 0) on invalid pixels.

    Returns
    -------
    shifts : ndarray
        Shift vector (in pixels) required to register ``target_image`` with
        ``src_image``.  Axis ordering is consistent with numpy (e.g. Z, Y, X)
	
	References
	----------
	.. [#] Dirk Padfield. "Masked Object Registration in the Fourier Domain". 
		IEEE Transactions on Image Processing, vol.21(5), pp. 2706-2718, 2012. 
	"""
	if mask is None:
		mask = np.ones_like(target_image)

	xcorr = masked_xcorr(target_image, src_image, m1 = mask, m2 = mask)

	maxima = np.unravel_index(np.argmax(np.abs(xcorr)), xcorr.shape)
	midpoints = np.array([np.fix(axis_size / 2) for axis_size in xcorr.shape])
	shifts = np.array(maxima, dtype = np.float)
	shifts[shifts > midpoints] -= np.array(xcorr.shape)[shifts > midpoints]

	return shifts

def masked_xcorr(arr1, arr2, m1, m2):
	""" 
	Normalized cross-correlation between two images with invalid pixels. 
	
	Parameters
	----------
	arr1, arr2 : `~numpy.ndarray`, shape (N,M)
		Images to be cross-correlated
	m1, m2 : `~numpy.ndarray`, shape (N,M)
		Masks of `arr1` and `arr2` respectively. The masks should evaluate to `False`
		(or 0) on invalid pixels.

	Returns
	-------
	out : `~numpy.ndarray`, shape
		Masked, normalized cross-correlation.

	References
	----------
	.. [#] Dirk Padfield. "Masked Object Registration in the Fourier Domain". 
		IEEE Transactions on Image Processing, vol.21(5), pp. 2706-2718, 2012. 
	"""
	m1, m2 = np.array(m1, dtype = np.float), np.array(m2, dtype = np.float)
	arr1, arr2 = np.asfarray(arr1), np.asfarray(arr2)

	# Contrary to the reference above, we do not compute the fft up to the 'full' size
	# In my experience this has not been a problem.
	fft2 = fftpack.fft2
	ifft2 = fftpack.ifft2
	
	F1 = fft2(arr1 * m1)
	F2s = fft2(arr2 * m2)
	np.conjugate(F2s, out = F2s)

	M1 = fft2(m1)
	M2s = fft2(m2)
	np.conjugate(M2s, out = M2s)

	# I have noticed no clear performance boost by storing
	# repeated calculation (e.g. ifft2(M1 * M2s)); however, the following
	# is already hard enough to read...
	numerator = ifft2(F1 * F2s)
	numerator -= ifft2(F1 * M2s) * ifft2(M1 * F2s) / ifft2(M1 * M2s)

	denominator = ifft2(fft2(arr1*arr1) * M2s) - (ifft2(F1 * M2s))**2/ifft2(M1 * M2s)
	denominator *= ifft2(M1*fft2(arr2*arr2).conjugate()) - ifft2(M1 * F2s)**2/(ifft2(M1 * M2s))
	np.sqrt(denominator, out = denominator)

	xcorr = np.real(numerator / denominator)
	return np.clip(xcorr, -1, 1)

def register_translation(src_image, target_image, src_mask = None, target_mask = None, **kwargs):
	"""
	Masked image translation registration by cross-correlation. In the case where
	no mask is provided, this function falls back to `skimage.feature.register_translation`.

	Parameters
	----------
	src_image : ndarray
		Reference image.
	target_image : ndarray
		Image to register.  Must be same dimensionality as ``src_image``.
	src_mask : `~numpy.ndarray`, shape (N,M) or None, optional
		Masks of `src_image`. The masks should evaluate to `False`
		(or 0) on invalid pixels. If `src_mask` is not None
		but `target_mask` is None, `src_mask` is used for both.
	target_mask : `~numpy.ndarray`, shape (N,M), or None, optional
		Masks of `target_image`. The masks should evaluate to `False`
		(or 0) on invalid pixels.

	Returns
	-------
	shifts : ndarray
		Shift vector (in pixels) required to register ``target_image`` with
		``src_image``.  Axis ordering is consistent with numpy (e.g. Z, Y, X)

	References
	----------
	.. [#] Dirk Padfield. "Masked Object Registration in the Fourier Domain". 
		IEEE Transactions on Image Processing, vol.21(5), pp. 2706-2718, 2012. 
	"""
	# TODO: upsampling

	if (src_mask is None) and (target_mask is None):
		return skimage_register_translation(src_image = src_image, target_image = target_image,
											upsample_factor = 1, space = 'real')[0]
	
	if (src_mask is not None) and (target_mask is None):
		cross_correlation = masked_xcorr(src_image, target_image, src_mask, src_mask)
	else:
		cross_correlation = masked_xcorr(src_image, target_image, src_mask, target_mask)
	

	maxima = np.unravel_index(np.argmax(np.abs(cross_correlation)), cross_correlation.shape)
	midpoints = np.array([np.fix(axis_size / 2) for axis_size in cross_correlation.shape])

	shifts = np.array(maxima, dtype = np.float)
	shifts[shifts > midpoints] -= np.array(cross_correlation.shape)[shifts > midpoints]
	#shifts[:] = shifts[::-1]	#

	# If its only one row or column the shift along that dimension has no
	# effect. We set to zero.
	for dim in range(cross_correlation.ndim):
		if cross_correlation.shape[dim] == 1:
			shifts[dim] = 0

	return shifts