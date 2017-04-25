"""
Parallelization utilities.
"""
# TODO: document this
# TODO: parallel reduce?

from collections.abc import Sized
from functools import partial
import multiprocessing as mp

def pmap(func, iterable, args = tuple(), kwargs = {}, processes = None):
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
		Positional arguments of `function`. `function` will be called
		as function(i, *args, **kwargs) for each i in `iterable`.
	kwargs : dictionary, optional
		Keyword arguments of `function`. `function` will be called
		as function(i, *args, **kwargs) for each i in `iterable`.
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
	# No point in spinning up a process pool for a single process
	if processes == 1:
		return map(partial(func, *args, **kwargs), iterable)

	if not isinstance(iterable, Sized):
		iterable = tuple(iterable)
	
	with mp.Pool(processes) as pool:
		# Best chunking is largest possible chunking
		chunksize = max(1, int(len(iterable)/pool._processes))
		
		map_func = pool.map
		if args:
			map_func = pool.starmap
			iterable = ((i,) + args for i in iterable)

		return map_func(func = partial(func, **kwargs), 
						iterable = iterable, 
						chunksize = chunksize)