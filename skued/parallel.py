"""
Parallelization utilities.

Functional programming-style `map` and `reduce` procedures are easily
parallelizable. The speed gain of parallelization can offset the
cost of spawning multiple processes for large iterables.
"""
# TODO: document this

from collections.abc import Sized
from functools import partial, reduce
import multiprocessing as mp

def chunked(iterable, chunksize = 1):
	"""
	Generator yielding multiple iterables of length 'chunksize'.

	Parameters
	----------
	iterable : iterable
		Must be a sized iterable, otherwise the iterable is consumed.
	chunksize : int, optional
	"""
	if not isinstance(iterable, Sized):
		iterable = tuple(iterable)
	length = len(iterable)
	for ndx in range(0, length, chunksize):
		yield iterable[ndx:min(ndx + chunksize, length)]

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
	added benefit of using `args` and `kwargs`.
	"""
	func = partial(func, *args, **kwargs)

	if processes == 1:
		return reduce(func, iterable)

	if not isinstance(iterable, Sized):
		iterable = tuple(iterable)

	with mp.Pool(processes) as pool:
		chunksize = max(1, int(len(iterable)/pool._processes))
		
		res = pool.imap_unordered(partial(reduce, func), tuple(chunked(iterable, chunksize)))
		return reduce(func, res)

def pmap(func, iterable, args = tuple(), kwargs = dict(), processes = None):
	"""
	Parallel application of a function with keyword arguments. Note
	that generators are consumed; hence, this implementation is
	not suitable for 'streaming' large objects through a pipeline...yet.

	Parameters
	----------
	func : callable
		Function to be applied to every element of `iterable`.
	iterable : iterable
		Iterable of items to be mapped. Generators are consumed.
	args : tuple
		Positional arguments of `function`.
	kwargs : dictionary, optional
		Keyword arguments of `function`.
	processes : int or None, optional
		Number of processes to use. If `None`, maximal number of processes
		is used. 

	Returns
	-------
	out : iterable
		Mapped values.

	Notes
	-----
	If `processes` is 1, `pmap` reduces to `map`, with the added benefit of
	of using `args` and `kwargs`
	"""
	if processes == 1:
		return map(partial(func, *args, **kwargs), iterable)

	if not isinstance(iterable, Sized):
		iterable = tuple(iterable)
	
	with mp.Pool(processes) as pool:
		chunksize = max(1, int(len(iterable)/pool._processes))

		map_func = pool.map
		if args:
			map_func = pool.starmap
			iterable = ((i,) + args for i in iterable)

		return map_func(func = partial(func, **kwargs), 
						iterable = iterable, 
						chunksize = chunksize)