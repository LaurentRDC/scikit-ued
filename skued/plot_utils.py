# -*- coding: utf-8 -*-
"""
Plotting utilities, including color manipulation. This is especially useful when
plotting time-series.
"""
from collections import Iterable
from colorsys import hsv_to_rgb

def linspace(start, stop, num):
	""" 
	Generate linear space 

    Parameters
    ----------
    start : float
        The starting value of the sequence.
    stop : float
        The end value of the sequence.
    num : int, optional
        Number of samples to generate.
	
	Yields
	------
	val : float
	"""
	step = (stop - start)/num

	val = start
	for _ in range(num):
		yield val
		val += step

def spectrum_colors(num_colors):
	"""
	Generates a set of RGB colors corresponding to the visible spectrum (i.e. the rainbow). 
	These colors can be used by Matplotlib, PyQtGraph, Qt, and other plotting/graphics libraries.
    
	Parameters
	----------
	num_colors : int or iterable
		Number of colors to generate. Alternatively, if `num_colors` is an iterable of numbers, 
		then the number of colors is deduced with specific spacing.
    
	Yields
	------
	color : (R,G,B) tuple.
		R, G, B values in a tuple, from 0.0 to 1.0.
	"""
	if isinstance(num_colors, int):
		num_colors = range(num_colors)
	num_colors = list(num_colors)

	if len(num_colors) == 1:
		hue_values = [0]
	else:
		mi = min(num_colors)
		ma = max(num_colors) - mi
		hue_values = reversed([(n - mi)/ma for n in num_colors])

	# Scale so that the maximum is 'purple' (hue = 0.8)
	# otherwise plots don't look as good
	yield from map(lambda hue: hsv_to_rgb(0.8*hue, s = 0.7, v = 0.9), hue_values)

def rgb_sweep(num_colors, source, dest):
	"""
	Generate a set of RGB colors as a linear sweep between two RGB colors `source` and `dest`.
	These colors can be used by Matplotlib, PyQtGraph, Qt, and other plotting/graphics libraries.

	Parameters
	----------
	num_colors : int
		Number of colors to generate.
	source : 3-tuple of floats
		Source color in RGB space. Values should be between 0.0 and 1.0.
	dest : 3-tuple of floats
		Destination color in RGB space. Values should be between 0.0 and 1.0.
	
	Yields
	------
	color : (R,G,B) tuple.
		R, G, B values in a tuple, from 0.0 to 1.0.
	"""	
	red_src, grn_src, blu_src = source
	red_dst, grn_dst, blu_dst = dest

	reds = linspace(red_src, red_dst, num = num_colors)
	greens = linspace(grn_src, grn_dst, num = num_colors)
	blues = linspace(blu_src, blu_dst, num = num_colors)

	yield from zip(reds, greens, blues)