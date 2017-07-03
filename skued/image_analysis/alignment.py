# -*- coding: utf-8 -*-
"""
Module concerned with alignment of diffraction images
"""
from itertools import product
import numpy as np
from skimage.feature import register_translation
from warnings import warn

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
	
	See also
	--------
	ialign
		generator of aligned images
	"""
	shifts, *_ = register_translation(target_image = _crop_to_half(image), 
									  src_image = _crop_to_half(reference), 
									  upsample_factor = 4, space = 'real')
	return shift_image(image, -shifts, fill_value = fill_value)

def _crop_to_half(image):
	nrows, ncols = np.array(image.shape)/4
	return image[int(nrows):-int(nrows), int(ncols):-int(ncols)]

def diff_register(image, reference, mask = None, search_space = 10):
	"""
	Register translation of diffraction patterns, using a 
	global optimization approach.
	
	Parameters
	----------
	image : iterable
		Iterable of ndarrays of shape (N,M)
	reference : `~numpy.ndarray`
		This is the reference image to which `image` will be aligned. 
	mask : `~numpy.ndarray` or None, optional
		Mask that evaluates to True on invalid pixels.
	search_space : int, optional
		Size of the domain (in pixels) over which a possible solution
		is computed. 
	
	Returns
	-------
	shift : `~numpy.ndarray`, shape (2,)
		Shift in rows and columns. Compatible with `shift_image`
	"""
	image = np.asfarray(image)
	reference = np.asfarray(reference)

	if mask is None:
		mask = np.zeros_like(image, dtype = np.bool)

	cropped = _crop_to_half(image)
	cropped_ref = _crop_to_half(reference)
	cropped_mask = _crop_to_half(mask)

	# I think that using nansum will be faster that alternatively
	# masking the shifted array everytime (which involves shifting the mask)
	cropped[cropped_mask] = np.nan
	cropped_ref[cropped_mask] = np.nan

	shifted = np.empty_like(cropped)	# preallocation for speed
	def cost(shift):
		""" Square of residual difference between im and reference 
		if im was shifted by i rows and j columns """
		shifted[:] = shift_image(cropped, shift, fill_value = np.nan)
		return np.nansum(np.square(shifted - cropped_ref))
	
	shift_space = list(product(range(-search_space, search_space + 1), range(-search_space, search_space + 1)))
	solutions = dict(zip(shift_space, map(cost, shift_space)))
	minimum = min(solutions.values())

	# Handling possible multiple minima is clunky af
	minimum_solutions = dict(filter(lambda item: item[-1] == minimum, solutions.items()))
	shift = sum(map(np.array, minimum_solutions.keys()))/len(minimum_solutions)
	return -1 * shift