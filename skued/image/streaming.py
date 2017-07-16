"""
Streaming operations on arrays/images
=====================================
"""

from functools import partial
from itertools import repeat
from math import sqrt

import numpy as np

from . import align

def ialign(images, reference = None, mask = None, fill_value = 0.0):
	"""
    Generator of aligned diffraction images.

    Parameters
    ----------
    images : iterable
        Iterable of ndarrays of shape (N,M)
    reference : `~numpy.ndarray`, shape (M,N)
        Images in `images` will be aligned onto the `reference` image. If
        'reference' is None (default), the first image in the 'images' stream
        is used as a reference
    mask : `~numpy.ndarray` or None, optional
        Mask that evaluates to True on invalid pixels.
    fill_value : float, optional
        Edges will be filled with `fill_value` after alignment.

    Yields
    ------
    aligned : `~numpy.ndarray`
        Aligned image. If `reference` is None, the first aligned image is the reference.

    See Also
    --------
    skued.image.align : align a single diffraction pattern onto a reference.
	"""
	images = iter(images)
	
	if reference is None:
		reference = next(images)
		yield reference

	yield from map(partial(align, reference = reference, mask = mask, fill_value =  fill_value), images)

def iaverage(images, weights = None):
    """ 
    Streaming average of diffraction images. This is equivalent
    to `numpy.average(axis = 2)` for a stack of images.

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
    numpy.average : (weighted) average for dense arrays
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
    Streaming variance of a set of images, per pixel. Weights are also supported.
    This is equivalent to calling `numpy.var(axis = 2)` on a stack of images.
    
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
    var: `~numpy.ndarray`
        Variance on a per-pixel basis. 
    
    See Also
    --------
    numpy.var : variance calculation for dense arrays. Weights are not supported.
    
    References
    ----------
    .. [#] D. H. D. West, Updating the mean and variance estimates: an improved method.
        Communications of the ACM Vol. 22, Issue 9, pp. 532 - 535 (1979)
    """
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
    Streaming standard deviation of images. Weights are also supported.
    This is equivalent to calling `numpy.std(axis = 2)` on a stack of images.

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
    std: `~numpy.ndarray`
        Standard deviation on a per-pixel basis.

    See Also
    --------
    numpy.std : standard deviation calculation of dense arrays. Weights are not supported.
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
    
    See Also
    --------
    scipy.stats.sem : standard error in the mean of dense arrays.
    """
    # TODO: include weights
    for k, std in enumerate(istd(images, ddof = ddof), start = 1):
        yield std / sqrt(k) 
