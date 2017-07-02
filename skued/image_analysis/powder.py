# -*- coding: utf-8 -*-
"""
Image manipulation of powder diffraction
"""
import numpy as np
from .alignment import diff_register

def powder_center(image, mask = None, search_space = 30):
	"""
	Finds the center of a powder diffraction pattern by comparing the
	correlation between the input and its image.

	Parameters  
	----------
	image : `~numpy.ndarray`, ndim 2
		Image of a powder pattern
	mask : `~numpy.ndarray` or None, optional
		Mask of `image`. The mask should evaluate to `True`
		(or 1) on invalid pixels. If None (default), no mask
		is used.
	search_space : int, optional
		Size of the domain (in pixels) over which a possible solution
		is computed. 

	Returns
	-------
	center : 2-tuple
		Center of the powder pattern. The center is returned in the format
		relevant for array manipulations (center = [row, column] instead of 
		center = [x,y]).
	"""
	if mask is None:
		mask = np.zeros_like(image, dtype = np.bool)

	shift = diff_register(image, np.rot90(image, k = 2), mask = mask * np.rot90(mask, k = 2), search_space = search_space)
	midpoints = np.array([int(axis_size / 2) for axis_size in image.shape])
	
	# I have found that there is always a residual (0.5, 0.5)
	center = shift[::-1]/2 + midpoints - np.array([0.5, 0.5])
	return tuple(center)

def angular_average(image, center, mask = None, extras = None, angular_bounds = None):
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
	angular_bounds : 2-tuple or None, optional
		If not None, the angles between first and second elements of `angular_bounds`
		(inclusively) will be used for the average. Angle bounds are specified in degrees.
		0 degrees is defined as the positive x-axis. Angle bounds outside [0, 360) are mapped back
		to [0, 360).

	Returns
	-------
	radius : `~numpy.ndarray`, shape (N,)
		Radius of the average [px]
	average : `~numpy.ndarray`, shape (N,)
		Angular-average of the array.
	"""
	#TODO: axis parameter? so batches of angular averages can be computed
	#	   at the same time.
	image = np.asarray(image, dtype = np.float)
	
	if mask is None:
		mask = np.zeros_like(image, dtype = np.bool)
	
	xc, yc = center  
	
	#Create meshgrid and compute radial positions of the data
	# TODO: is there no way to use rint and get a dtype np.int at the end?
	#		astype() takes about 20% of computing time of this function.
	Y, X = np.ogrid[0:image.shape[0],0:image.shape[1]]
	R = np.rint(np.sqrt( (X - xc)**2 + (Y - yc)**2 )).astype(np.int)

	if angular_bounds:
		mi, ma = angular_bounds
		mi, ma = mi % 360, ma % 360
		angles = np.rad2deg(np.arctan2(Y - yc, X - xc))
		mask[np.logical_not(np.logical_and(mi <= angles, angles <= ma))] = True

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
	
	return np.arange(0, radial_intensity.size), radial_intensity