# -*- coding: utf-8 -*-
"""
Image manipulation involving symmetry
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import fftconvolve, find_peaks_cwt
from skimage.filters import gaussian, threshold_local

import matplotlib.pyplot as plt

_PC_CACHE = dict()
def powder_center(image):
	"""
	Finds the center of a powder diffraction pattern by comparing the
	correlation between the input and its image.

	Parameters  
	----------
	image : `~numpy.ndarray`, ndim 2
		Image of a powder pattern

	Returns
	-------
	center : 2-tuple
		Center of the powder pattern. The center is returned in the format
		relevant for array manipulations (center = [row, column] instead of 
		center = [x,y]).
	"""
	# The correlation 'background' only depends on the image
	# shape; therefore, caching provides an important speedup
	if image.shape not in _PC_CACHE:
		mask = np.ones_like(image)
		_PC_CACHE[image.shape] = fftconvolve(mask, mask)

	# Correlation of image and its mirror is a convolution
	# TODO: explicitly compute convolution, since only one FFT
	# 		has to be calculated.
	corr = fftconvolve(image, image)
	corr /= _PC_CACHE[image.shape]

	# ignore edges because of artifacts from fftconvolve
	edge_size = int(min(image.shape) / 10)
	corr[:edge_size,:] = 0
	corr[-edge_size:, :] = 0
	corr[:, -edge_size:] = 0
	corr[:, :edge_size] = 0

	# Noise in the image will obfuscate the symmetry peak
	# Therefore, we have to identify where the (sharp) peak is
	thresh = threshold_local(corr, block_size = 101)
	corr[corr <= thresh] = 0

	full_center = np.array(np.unravel_index(np.argmax(corr), dims = corr.shape))
	return tuple(full_center/2)

def angular_average(image, center, mask = None, extras = None):
	"""
	This function returns an angularly-averaged pattern computed from a diffraction pattern, 
	e.g. polycrystalline diffraction.
    
	Parameters
	----------
	arr : array_like, shape (M, N)
		Array or image.
	center : array_like, shape (2,)
		coordinates of the center (in pixels).
	mask : `~numpy.ndarray` or None, optional
		Evaluates to True on invalid elements of array.
	extras : dict-like or None, optional
		if not None, this dict-like object will be updated with: 

			extras['error'] : `~numpy.ndarray`
				standard error in mean across radii.

	Returns
	-------
	radius : `~numpy.ndarray`, shape (N,)
		Radius of the average [px]
	average : `~numpy.ndarray`, shape (N,)
		Angular-average of the array.
	"""
	#TODO: N dimensions. The only problem is the generation of an arbitrary
	#	   meshgrid on the image.
	#TODO: axis parameter? so batches of angular averages can be computed
	#	   at the same time.
	image = np.asarray(image, dtype = np.float)
	
	if mask is None:
		mask = np.zeros_like(image, dtype = np.bool)
	
	xc, yc = center  
	
	#Create meshgrid and compute radial positions of the data
	X, Y = np.meshgrid(np.arange(0, image.shape[0]), 
					   np.arange(0, image.shape[1]))
	R = np.rint(np.sqrt( (X - xc)**2 + (Y - yc)**2 )).astype(np.int)

	valid = np.logical_not(mask)
	R_v, image_v = R[valid].ravel(), image[valid].ravel()
	px_bin = np.bincount(R_v, weights = image_v)
	r_bin = np.bincount(R_v)

	# np.bincount will start counting at 0. We ignore the leading zeroes
	nz = r_bin > 0.0
	radial_intensity = px_bin[nz]/r_bin[nz]

	# Update the extras dictionary if provided:
	# Error as the standard error in the mean, at each pixel
	# Standard error = std / sqrt(N)
	# std = sqrt(var - mean**2)
	if extras is not None:
		var_bin = np.bincount(R_v, weights = image_v**2)[nz]/r_bin[nz]
		radial_intensity_error = np.sqrt(var_bin - radial_intensity**2)/np.sqrt(r_bin[nz])
		extras.update({'error':radial_intensity_error})
	
	return np.unique(R_v), radial_intensity