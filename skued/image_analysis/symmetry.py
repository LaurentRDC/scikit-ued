"""
Image manipulation involving symmetry
"""

import numpy as np

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

		``extras['error'] : `~numpy.ndarray` ``
			Standard error in mean across radii.

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
	
	xc, yc = center     #Center coordinates  
	
	#Create meshgrid and compute radial positions of the data
	X, Y = np.meshgrid(np.arange(0, image.shape[0]), 
					   np.arange(0, image.shape[1]))
	R = np.rint(np.sqrt( (X - xc)**2 + (Y - yc)**2 )).astype(np.int)

	# Replace all values in the image corresponding to irrelevant# data by 0: 
	# this way, it will never count in any calculation because image
	# values are used as weights in numpy.bincount
	image[mask] = 0
	
	# Angular average
	px_bin = np.bincount(R.ravel(), weights = image.ravel())
	r_bin = np.bincount(R.ravel())
	radial_intensity = px_bin/r_bin

	# Update the extras dictionary if provided:
	# Error as the standard error in the mean, at each pixel
	# Standard error = std / sqrt(N)
	# std = sqrt(var - mean**2)
	if extras is not None:
		var_bin = np.bincount(R.ravel(), weights = image.ravel()**2)/r_bin
		radial_intensity_error = np.sqrt(var_bin - radial_intensity**2)/np.sqrt(r_bin)
		extras.update({'error':radial_intensity_error})

	return np.unique(R), radial_intensity