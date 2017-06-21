"""
Array utility functions
"""

from itertools import repeat
from numpy import concatenate, swapaxes

def repeated_array(arr, num, axes = -1):
    """
    Create a composite array from repeated copies of an array
    
    Parameters
    ----------
    arr : ndarray
    
    num : int or iterable of ints
        Number of copies per axis. If provided as tuple, must be the same length
        as 'axes' parameter. In case of `num` being 0 or an empty iterable, 
        the inpur `arr` is returned.
    axes : int or iterable of ints
        Axis/axes over which to copy.
    
    Returns
    -------
    out : ndarray
    
    Raises
    ------
    ValueError
        If num and axes are tuples of different lengths.
    """
    if not num:
        return arr
    
    if isinstance(num, int): num = (num,) 
    if isinstance(axes, int): axes = (axes,)
    
    if len(num) != len(axes):
        raise ValueError('num and axes must have the same length')
    
    composite = concatenate(tuple(repeat(arr, times = num[0])), axis = axes[0])

    if len(num) > 1:
        for n, ax in zip(num[1:], axes[1:]):
            composite = concatenate(tuple(repeat(composite, times = n)), axis = ax)
    
    return composite

def mirror(arr, axes = None):
    """ 
    Reverse array over many axes. Generalization of arr[::-1] for many dimensions.

    Parameters
    ----------
    arr : `~numpy.ndarray`
        Array to be reversed
    axes : int or tuple or None, optional
        Axes to be reversed. Default is to reverse all axes.
    
    Returns
    -------
    out : 
    """
    if axes is None:
        axes = range(arr.ndim)

    elif isinstance(axes, int):
        axes = (axes,)
    
    # arr[::-1] reverses in the first axis direction
    for axis in axes:
        arr = swapaxes(swapaxes(arr, 0, axis)[::-1], 0, axis)
    
    return arr