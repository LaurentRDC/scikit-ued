# -*- coding: utf-8 -*-
"""
Parallelization utilities
=========================

Functional programming-style `map` and `reduce` procedures are easily
parallelizable. The speed gain of parallelization can offset the
cost of spawning multiple processes for large iterables.
"""
import multiprocessing as mp
from collections.abc import Sized
from functools import partial, reduce

from .iter_utils import chunked

def preduce(func, iterable, args = tuple(), kwargs = dict(), processes = None):
	"""
	Parallel application of the reduce function, with keyword arguments.

	Parameters
	----------
	func : callable
		Function to be applied to every element of `iterable`.
	iterable : iterable
		Iterable of items to be reduced. Generators are consumed.
	args : tuple
		Positional arguments of `function`.
	kwargs : dictionary, optional
		Keyword arguments of `function`.
	processes : int or None, optional
		Number of processes to use. If `None`, maximal number of processes
		is used. 

	Returns
	-------
	reduced : object

	Notes
	-----
	If `processes` is 1, `preduce` is equivalent to functools.reduce with the
	added benefit of using `args` and `kwargs`, but `initializer` is not supported.
	"""
	func = partial(func, *args, **kwargs)

	if processes == 1:
		return reduce(func, iterable)

	with mp.Pool(processes) as pool:
		if isinstance(iterable, Sized):
			chunksize = max(1, int(len(iterable)/pool._processes))
		else:
			chunksize = 1
		
		res = pool.imap_unordered(partial(reduce, func), tuple(chunked(iterable, chunksize)))
		return reduce(func, res)

def pmap(func, iterable, args = tuple(), kwargs = dict(), processes = None):
	"""
	Parallel application of a function with keyword arguments. Based on 
	multiprocessing.Pool.imap.

	Parameters
	----------
	func : callable
		Function to be applied to every element of `iterable`.
	iterable : iterable
		Iterable of items to be mapped.
	args : tuple, optional
		Positional arguments of `function`.
	kwargs : dictionary, optional
		Keyword arguments of `function`.
	processes : int or None, optional
		Number of processes to use. If `None`, maximal number of processes
		is used. 

	Yields
	------
	Mapped values.

	Notes
	-----
	If `processes` is 1, `pmap` reduces to `map`, with the added benefit of
	of using `kwargs`
	"""
	func = partial(func, *args, **kwargs)

	if processes == 1:
		yield from map(func, iterable)
	
	with mp.Pool(processes) as pool:
		if isinstance(iterable, Sized):
			chunksize = max(1, int(len(iterable)/pool._processes))
		else:
			chunksize = 1

		yield from pool.imap(func = func, iterable = iterable, chunksize = chunksize)
