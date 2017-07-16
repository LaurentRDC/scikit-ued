"""
Streaming operations on arrays/images
=====================================
"""
from collections import deque
from functools import partial
from itertools import repeat
from math import sqrt

import numpy as np

from . import align


def ialign(images, reference = None, fill_value = 0.0):
	"""
	Generator of aligned diffraction images.

	Parameters
	----------
	images : iterable
		Iterable of ndarrays of shape (N,M)
	reference : `~numpy.ndarray` or None, optional
		If not None, this is the reference image to which all images will be aligned. Otherwise,
		images will be aligned to the first element of the iterable 'images'. 
	fill_value : float, optional
		Edges will be filled with `fill_value` after shifting.
    
	Yields
	------
	aligned : ndarray, ndim 2
		Aligned image

	Notes
	-----
	Diffraction images exhibit high symmetry in most cases, therefore images
	are cropped to a quarter of their size before alignment.
	"""
	images = iter(images)
	
	if reference is None:
		reference = next(images)
		yield reference

	yield from map(partial(align, reference = reference), images)

def iaverage(images, weights = None):
    """ 
    Streaming average of diffraction images. This generator can be used to 
    observe a live averaging.

    Parameters
    ----------
    images : iterable of ndarrays
        Images to be averaged. This iterable can also a generator.
    weights : iterable of ndarray, iterable of floats, or None, optional
        Iterable of weights associated with the values in each item of `images`. 
        Each value in an element of `images` contributes to the average 
        according to its associated weight. The weights array can either be a float
        or an array of the same shape as any element of `images`. If weights=None, 
        then all data in each element of `images` are assumed to have a weight equal to one.
    
    Yields
    ------
    avg: `~numpy.ndarray`
        Weighted average. 
    
    See Also
    --------
    numpy.average : average for dense arrays
    """
    images = iter(images)
    
    if weights is None:
        weights = repeat(1.0)
    weights = iter(weights)

    sum_of_weights = np.array(next(weights), copy = True)
    weighted_sum = np.array(next(images) * sum_of_weights, copy = True)
    yield weighted_sum/sum_of_weights

    for image, weight in zip(images, weights):

        sum_of_weights += weight
        weighted_sum += weight * image
        #print(sum_of_weights)
        yield weighted_sum/sum_of_weights

def ivar(images, ddof = 1, weights = None):
    """ 
    Streaming variance of a set of images, per pixel. This is equivalent to
    calling `numpy.var` with `ddof = 1`.

    Parameters
    ----------
    images : iterable of ndarrays
        Images to be combined. This iterable can also a generator.
    ddof : int, optional
        Means Delta Degrees of Freedom.  The divisor used in calculations
        is ``N - ddof``, where ``N`` represents the number of elements.
        By default `ddof` is one.
    weights : iterable of ndarray, iterable of floats, or None, optional
        Iterable of weights associated with the values in each item of `images`. 
        Each value in an element of `images` contributes to the variance 
        according to its associated weight. The weights array can either be a float
        or an array of the same shape as any element of `images`. If weights=None, 
        then all data in each element of `images` are assumed to have a weight equal to one.
    
    Yields
    ------
    sem: `~numpy.ndarray`
        Variance on a per-pixel basis. 
    
    References
    ----------
    .. [#] D. Knuth, The Art of Computer Programming 3rd Edition, Vol. 2, p. 232
    """
    # TODO: weighted online variance 
    # See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
    images = iter(images)

    if weights is None:
        weights = repeat(1.0)
    weights = iter(weights)
    sum_of_weights = np.array(next(weights), copy = True)

    first = next(images)
    old_mean = new_mean = np.array(first, copy = True)
    old_S = new_S = np.zeros_like(first, dtype = np.float)
    yield np.zeros_like(first)  # No error if no averaging
    
    for image, weight in zip(images, weights):

        sum_of_weights += weight

        _sub = weight * (image - old_mean)
        new_mean[:] = old_mean + _sub/sum_of_weights
        new_S[:] = old_S + _sub*(image - new_mean)

        # In case there hasn't been enough measurements yet,
        # yield zeros.
        if np.any(sum_of_weights - ddof <= 0):
            yield np.zeros_like(first)
        else:
            yield new_S/(sum_of_weights - ddof) # variance = S / k-1, sem = std / sqrt(k)    

        old_mean[:] = new_mean
        old_S[:] = new_S

def istd(images, ddof = 1, weights = None):
    """ 
    Streaming standard deviation of images. This is equivalent to
    calling `numpy.std(axis = 2)` on a stack of images.

    Parameters
    ----------
    images : iterable of ndarrays
        Images to be combined. This iterable can also a generator.
    ddof : int, optional
        Means Delta Degrees of Freedom.  The divisor used in calculations
        is ``N - ddof``, where ``N`` represents the number of elements.
        By default `ddof` is one.
    weights : iterable of ndarray, iterable of floats, or None, optional
        Iterable of weights associated with the values in each item of `images`. 
        Each value in an element of `images` contributes to the standard deviation
        according to its associated weight. The weights array can either be a float
        or an array of the same shape as any element of `images`. If weights=None, 
        then all data in each element of `images` are assumed to have a weight equal to one.
    
    Yields
    ------
    sem: `~numpy.ndarray`
        Standard deviation on a per-pixel basis.
    """
    yield from map(np.sqrt, ivar(images, ddof = ddof, weights = weights))

def isem(images, ddof = 1):
    """ 
    Streaming standard error in the mean (SEM) of images. This is equivalent to
    calling `scipy.stats.sem(axis = 2)` on a stack of images.

    Parameters
    ----------
    images : iterable of ndarrays
        Images to be combined. This iterable can also a generator.
    ddof : int, optional
        Means Delta Degrees of Freedom.  The divisor used in calculations
        is ``N - ddof``, where ``N`` represents the number of elements.
        By default `ddof` is one.
    
    Yields
    ------
    sem: `~numpy.ndarray`
        Standard error in the mean. 
    
    See also
    --------
    scipy.stats.sem : standard error in the mean of dense arrays.
    """
    for k, std in enumerate(istd(images, ddof = ddof), start = 1):
        yield std / sqrt(k) 
