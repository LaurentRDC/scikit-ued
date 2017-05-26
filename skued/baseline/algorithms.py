# -*- coding: utf-8 -*-
"""
Algorithms based on the dual-tree complex wavelet transform.

Author : Laurent P. Rene de Cotret
"""
from collections.abc import Iterable
import pywt
from . import dtcwt, idtcwt
import numpy as np
from functools import partial
from warnings import warn

def _iterative_baseline(array, max_iter, mask, background_regions, axes, approx_rec_func, func_kwargs):
	""" Base function for iterative baseline determination """
	array = np.asarray(array, dtype = np.float)

	if isinstance(axes, int):
		axes = (axes,)

	approx_rec = partial(approx_rec_func, **func_kwargs)

	# Since most wavelet transforms only works on even-length signals, we might have to extend.
	# See numpy.pad() docs for a formatting of the padding tuple constructed below
	original_shape = array.shape
	padding = [(0,0) for dim in range(array.ndim)]     # e.g. 2D : padding = [ (0,0), (0,0) ]
	for axis in axes:
		if original_shape[axis] % 2 == 1:
			padding[axis] = (0, 1) 
	array = np.pad(array, tuple(padding), mode = 'edge')

	if mask is None:
		mask = np.zeros_like(array, dtype = np.bool)
	
	signal = np.array(array, copy = True)
	background = np.zeros_like(array, dtype = np.float)
	for i in range(max_iter):
		
		# Make sure the background values are equal to the original signal values in the
		# background regions
		for index in background_regions:
			signal[index] = array[index]

		# Wavelet reconstruction using approximation coefficients
		background[:] = approx_rec(signal)
		
		# Modify the signal so it cannot be more than the background
		# This reduces the influence of the peaks in the wavelet decomposition
		signal[signal > background] = background[signal > background]
	
	# The background should be identically 0 where the data points are invalid
	background[mask] = 0 

	# Readjust size for odd input signals
	# If axis was padded, trim
	if background.shape != original_shape:
		for axis in axes:
			background = np.swapaxes(np.swapaxes(background, 0, axis)[:-1], 0, axis)
	return background

def _dt_approx_rec(array, first_stage, wavelet, level, axis = -1):
    """
    Approximate reconstruction of a signal/image using the dual-tree approach.
    
    Parameters
    ----------
    array : array-like
        Array to be decomposed.
    level : int or 'max'
        Decomposition level. A higher level will result in a coarser approximation of
        the input array. If the level is higher than the maximum possible decomposition level,
        the maximum level is used. If 'max' (default), the maximum possible decomposition level is used.
    first_stage : str, optional
        First-stage wavelet to use.
    wavelet : str, optional
        Complex wavelet to use in late stages.
    axis : int, optional
        Axis over which to compute the transform. Default is -1.
            
    Returns
    -------
    reconstructed : ndarray
        Approximated reconstruction of the input array.
    """
    coeffs = dtcwt(data = array, first_stage = first_stage, wavelet = wavelet, level = level, mode = 'constant', axis = axis)
    app_coeffs, *det_coeffs = coeffs
    
    det_coeffs = [np.zeros_like(det, dtype = np.complex) for det in det_coeffs]
    reconstructed = idtcwt(coeffs = [app_coeffs] + det_coeffs, first_stage = first_stage, wavelet = wavelet, mode = 'constant', axis = axis)
    return reconstructed

def _dwt_approx_rec(array, level, wavelet, axis):
	"""
	Approximate reconstruction of a signal/image. Uses the multi-level discrete wavelet 
	transform to decompose a signal or an image, and reconstruct it using approximate 
	coefficients only.
    
	Parameters
	----------
	array : array-like
		Array to be decomposed. Currently, only 1D and 2D arrays are supported.
		Only even-lengths signals long the axis.
	level : int or 'max' or None (deprecated)
		Decomposition level. A higher level will result in a coarser approximation of
		the input array. If the level is higher than the maximum possible decomposition level,
		the maximum level is used.
		If None, the maximum possible decomposition level is used.
	wavelet : str or Wavelet object
		Can be any argument accepted by PyWavelet.Wavelet, e.g. 'db10'
    axis : int, optional
        Axis over which to compute the transform. Default is -1.
            
	Returns
	-------
	reconstructed : ndarray
		Approximated reconstruction of the input array.
    
	Raises
	------    
	ValueError
		If input array has dimension > 2 
	"""
	if isinstance(axis, Iterable):
		axis = axis[0]

	array = np.asarray(array, dtype = np.float)
	if array.shape[axis] % 2 != 0:
		raise ValueError('Only even-length signals are supported')

	# Build Wavelet object
	if isinstance(wavelet, str):
		wavelet = pywt.Wavelet(wavelet)
	
	# Check maximum decomposition level
	# For 2D array, check the condition with shortest dimension min(array.shape). This is how
	# it is done in PyWavelet.wavedec2.
	max_level = pywt.dwt_max_level(data_len = array.shape[axis], filter_len = wavelet.dec_len)
	if level is None or level is 'max':
		level = max_level
	elif max_level < level:
		warn('Decomposition level {} higher than maximum {}. Maximum is used.'.format(level, max_level))
		level = max_level
	
	# By now, we are sure that the decomposition level will be supported.
	# Decompose the signal using the multilevel discrete wavelet transform
	coeffs = pywt.wavedec(data = array, wavelet = wavelet, level = level, mode = 'constant', axis = axis)
	app_coeffs, det_coeffs = coeffs[0], coeffs[1:]
	
	# Replace detail coefficients by 0; keep the correct length so that the
	# reconstructed signal has the same size as the (possibly upsampled) signal
	# The structure of coefficients depends on the dimensionality	
	# Reconstruct signal
	reconstructed = pywt.waverec([app_coeffs] + [None,]*len(det_coeffs), wavelet = wavelet, mode = 'constant', axis = axis)

	# Sometimes pywt.waverec returns a signal that is longer than the original signal
	if reconstructed.shape[axis] > array.shape[axis]:
		reconstructed = np.swapaxes(np.swapaxes(reconstructed, 0, axis)[:array.shape[axis]], 0, axis)
	return  reconstructed

def _dwt_approx_rec2(array, level, wavelet, axis):
	"""
	Approximate reconstruction of a signal/image. Uses the multi-level discrete wavelet 
	transform to decompose a signal or an image, and reconstruct it using approximate 
	coefficients only.
    
	Parameters
	----------
	array : array-like
		Array to be decomposed. Currently, only 1D and 2D arrays are supported.
		Only even-lengths signals long the axis.
	level : int or 'max' or None (deprecated)
		Decomposition level. A higher level will result in a coarser approximation of
		the input array. If the level is higher than the maximum possible decomposition level,
		the maximum level is used.
		If None, the maximum possible decomposition level is used.
	wavelet : str or Wavelet object
		Can be any argument accepted by PyWavelet.Wavelet, e.g. 'db10'
	axis : 2-tuple of ints
            
	Returns
	-------
	reconstructed : ndarray
		Approximated reconstruction of the input array.
    
	Raises
	------    
	ValueError
		If input array has dimension > 2 
	"""
	array = np.asarray(array, dtype = np.float)
	for ax in axis:
		if array.shape[ax] % 2 != 0:
			raise ValueError('Only even-length signals are supported')

	# Build Wavelet object
	if isinstance(wavelet, str):
		wavelet = pywt.Wavelet(wavelet)
	
	# Check maximum decomposition level
	# For 2D array, check the condition with shortest dimension min(array.shape). This is how
	# it is done in PyWavelet.wavedec2.
	max_level = pywt.dwt_max_level(data_len = min(array.shape[ax] for ax in axis), filter_len = wavelet.dec_len)
	if level is None:
		level = max_level
	elif max_level < level:
		warn('Decomposition level {} higher than maximum {}. Maximum is used.'.format(level, max_level))
		level = max_level
	
	# By now, we are sure that the decomposition level will be supported.
	# Decompose the signal using the multilevel discrete wavelet transform
	coeffs = pywt.wavedec2(data = array, wavelet = wavelet, level = level, mode = 'constant', axes = axis)
	app_coeffs, det_coeffs = coeffs[0], coeffs[1:]
	
	# Replace detail coefficients by 0; keep the correct length so that the
	# reconstructed signal has the same size as the (possibly upsampled) signal
	# The structure of coefficients depends on the dimensionality	
	# Reconstruct signal
	reconstructed = pywt.waverec2([app_coeffs] + [(None, None, None)]*len(det_coeffs), 
								  wavelet = wavelet, mode = 'constant', axes = axis)

	# Sometimes pywt.waverec returns a signal that is longer than the original signal
	for ax in axis:
		if reconstructed.shape[ax] > array.shape[ax]:
			reconstructed = np.swapaxes(np.swapaxes(reconstructed, 0, ax)[:-1], 0, ax)
	return reconstructed

def baseline_dt(array, max_iter, level = None, first_stage = 'sym6', wavelet = 'qshift1', 
				background_regions = [], mask = None, axis = -1):
	"""
	Iterative method of baseline determination based on the dual-tree complex wavelet transform.
	This function only works in 1D, along an axis. For baseline of 2D arrays, see baseline_dwt.
    
	Parameters
	----------
	array : ndarray, shape (M,N)
		Data with background. Can be either 1D signal or 2D array.
	max_iter : int
		Number of iterations to perform.
	level : int or 'max', optional
		Decomposition level. A higher level will result in a coarser approximation of
		the input signal (read: a lower frequency baseline). If 'max' (default), the maximum level
		possible is used.
	first_stage : str, optional
		Wavelet to use for the first stage. See skued.baseline.ALL_FIRST_STAGE for a list of suitable arguments
	wavelet : str, optional
		Wavelet to use in stages > 1. Must be appropriate for the dual-tree complex wavelet transform.
		See skued.baseline.ALL_COMPLEX_WAV for possible values.
	background_regions : iterable, optional
		Indices of the array values that are known to be purely background. Depending
		on the dimensions of array, the format is different:
        
		``array.ndim == 1``
			`background_regions` is a list of ints (indices) or slices
			E.g. >>> background_regions = [0, 7, 122, slice(534, 1000)]
          
		``array.ndim == 2``
			`background_regions` is a list of tuples of ints (indices) or tuples of slices
			E.g. >>> background_regions = [(14, 19), (42, 99), (slice(59, 82), slice(81,23))]
         
		Default is empty list.
    
	mask : ndarray, dtype bool, optional
		Mask array that evaluates to True for pixels that are invalid. 
	axis : int, optional
		Axis over which to compute the wavelet transform. Default is -1
    
	Returns
	-------
	baseline : ndarray, shape (M,N)
		Baseline of the input array.
        
	References
	----------
	.. [#] L. P. René de Cotret and B. J. Siwick, A general method for baseline-removal in ultrafast 
		   electron powder diffraction data using the dual-tree complex wavelet transform, Struct. Dyn. 4 (2016)
	"""
	return _iterative_baseline(array = array, max_iter = max_iter, background_regions = background_regions,
									mask = mask, axes = (axis,), approx_rec_func = _dt_approx_rec,
									func_kwargs = {'level': level, 'wavelet': wavelet,
													'first_stage': first_stage, 'axis': axis})

def baseline_dwt(array, max_iter, level = None, wavelet = 'sym6', background_regions = [], mask = None, axis = -1):
	"""
	Iterative method of baseline determination, based on the discrete wavelet transform. 
    
	Parameters
	----------
	array : ndarray
		Data with background.
	max_iter : int
		Number of iterations to perform.
	level : int or None, optional
		Decomposition level. A higher `level` will result in a coarser approximation of
		the input signal (read: a lower frequency baseline). If None (default), the maximum level
		possible is used.
	wavelet : PyWavelet.Wavelet object or str, optional
		Wavelet with which to perform the algorithm. See PyWavelet documentation
		for available values. Default is 'sym6'.
	background_regions : list, optional
		Indices of the array values that are known to be purely background. Depending
		on the dimensions of array, the format is different:
        
		``array.ndim == 1``
			`background_regions` is a list of ints (indices) or slices
			E.g. >>> background_regions = [0, 7, 122, slice(534, 1000)]
          
		``array.ndim == 2``
			`background_regions` is a list of tuples of ints (indices) or tuples of slices
			E.g. >>> background_regions = [(14, 19), (42, 99), (slice(59, 82), slice(81,23))]
         
		Default is empty list.
    
	mask : ndarray, dtype bool, optional
		Mask array that evaluates to True for pixels that are invalid. Useful to determine which pixels are masked
		by a beam block.
	axis : int or tuple, optional
		Axis over which to compute the wavelet transform. Can also be a 2-tuple of ints for 2D baseline
    
	Returns
	-------
	baseline : ndarray, shape (M,N)
		Baseline of the input array.
    
	References
	----------
	.. [#] Galloway et al. 'An Iterative Algorithm for Background Removal in Spectroscopy by Wavelet 
		   Transforms', Applied Spectroscopy pp. 1370 - 1376, September 2009.
	"""
	if isinstance(axis, int):
		axis = (axis,)
	
	approx_rec_func = {1: _dwt_approx_rec, 2: _dwt_approx_rec2}

	return _iterative_baseline(array, max_iter = max_iter, background_regions = background_regions, 
							   mask = mask, axes = axis, approx_rec_func = approx_rec_func[len(axis)], 
							   func_kwargs = {'level': level, 'wavelet': wavelet, 'axis': axis})