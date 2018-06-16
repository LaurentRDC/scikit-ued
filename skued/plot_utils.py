# -*- coding: utf-8 -*-
"""
Plotting utilities 
==================
"""
from colorsys import hsv_to_rgb

from npstreams import linspace, multilinspace


def _hex_to_rgb(value):
	""" Return an RGB tuple (float 0 to 1) from an hexadecimal 
	string representation. """
	value = value.lstrip('#')
	char_per_color = len(value) // 3

	red = int(value[0:char_per_color], base = 16)
	green = int(value[char_per_color:2*char_per_color], base = 16)
	blue = int(value[2*char_per_color::], base = 16) 

	# We want color values between 0.0 and 1.0
	max_bits = 16**char_per_color - 1
	return (red/max_bits, green/max_bits, blue/max_bits)

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
	source : 3-tuple of floats or str
		Source color in RGB space. Values should be between 0.0 and 1.0.
		if `source` is a string, it is assumed to be in hexadecimal representation.
	dest : 3-tuple of floats or str
		Destination color in RGB space. RGB tuple Values should be between 0.0 and 1.0.
		If `dest` is a string, it is assumed to be in hexadecimal representation.
	
	Yields
	------
	color : (R,G,B) tuple.
		R, G, B values in a tuple, from 0.0 to 1.0.
	"""	
	if isinstance(source, str):
		source = _hex_to_rgb(source)
	
	if isinstance(dest, str):
		dest = _hex_to_rgb(dest)

	yield from multilinspace(source, dest, num = num_colors)
