"""
Module concerned with alignment of diffraction images
"""

import numpy as np
from skimage.feature import register_translation

try:
	from numpy.fft_intel import fft2, ifft2
except ImportError:
	from scipy.fftpack import fft2, ifft2

non = lambda s: s if s < 0 else None
mom = lambda s: max(0, s)
def shift_image(arr, shift, fill_value = 0):
	""" 
	Shift an image at a 1-pixel resolution. Shift

	Parameters
	----------
	arr : ndarray
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

def align(images, reference = None, fill_value = 0.0):
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
    
	Yields
	------
	aligned : ndarray, ndim 2
		Aligned image

	Notes
	-----
	Diffraction images exhibit high symmetry in most cases, therefore images
	are cropped to a quarter of their size before alignment.
	"""
	images = iter(images)
	
	if reference is None:
		reference = next(images)
		yield reference

	cropped_ref = _crop_to_half(reference)

	def _align_im(image):
		shifts, *_ = register_translation(target_image = _crop_to_half(image), 
										  src_image = cropped_ref, 
										  upsample_factor = 4, 
										  space = 'real')
		return shift_image(image, -shifts, fill_value = fill_value)

	yield from map(_align_im, images)

def _crop_to_half(image):
	nrows, ncols = np.array(image.shape)/4
	return image[int(nrows):-int(nrows), int(ncols):-int(ncols)]