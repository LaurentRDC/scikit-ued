# -*- coding: utf-8 -*-
"""
Iterator/Generator utilities 
============================
"""
from collections import deque, Sized
from itertools import islice, count

def chunked(iterable, chunksize = 1):
	"""
	Generator yielding multiple iterables of length 'chunksize'.

	Parameters
	----------
	iterable : iterable
		Iterable to be chunked. 
	chunksize : int, optional
		Chunk size. 

	Yields
	------
	chunk : iterable
		Iterable of size `chunksize`. In special case of iterable not being
		divisible by `chunksize`, the last `chunk` might be smaller.
	"""
	# This looks ridiculously simple now, 
	# but I didn't always know about itertools
	iterable = iter(iterable)

	next_chunk = tuple(islice(iterable, chunksize))
	while next_chunk:	
		yield next_chunk
		next_chunk = tuple(islice(iterable, chunksize))
		

def linspace(start, stop, num, endpoint = True):
	""" 
	Generate linear space. This is sometimes more appropriate than
	using `range`.

	Parameters
	----------
	start : float
		The starting value of the sequence.
	stop : float
		The end value of the sequence.
	num : int, optional
		Number of samples to generate.
	endpoint : bool, optional
		If True (default), the endpoint is included in the linear space.

	Yields
	------
	val : float

	See also
	--------
	numpy.linspace : generate linear space as a dense array.
	"""
	# If endpoint are to be counted in,
	# step does not count the last yield
	if endpoint:
		num -= 1

	step = (stop - start)/num

	val = start
	for _ in range(num):
		yield val
		val += step

	if endpoint:
		yield stop

def multilinspace(start, stop, num, endpoint = True):
	""" 
	Generate multilinear space, for joining the values in two iterables.

	Parameters
	----------
	start : iterable
		The starting value. 
	stop : iterable
		The end value.
	num : int, optional
		Number of samples to generate.
	endpoint : bool, optional
		If True (default), the endpoint is included in the linear space.

	Yields
	------
	val : tuple
		Tuple of the same length as start and stop

	See also
	--------
	linspace : generate a linear space between two numbers
	"""	
	start, stop = tuple(start), tuple(stop)
	if len(start) != len(stop):
		raise ValueError('start and stop must have the same length')

	spaces = tuple(linspace(a, b, num = num, endpoint = endpoint) for a, b in zip(start, stop))
	yield from zip(*spaces)

def last(stream):
	""" Retrieve the last item from a stream/iterator. Generators are consumed. 
	If empty stream, returns None """
	# Wonderful idea from itertools recipes
	# https://docs.python.org/3.6/library/itertools.html#itertools-recipes
	try:
		return deque(stream, maxlen = 1)[0]
	except IndexError:	# Empty stream
		return None