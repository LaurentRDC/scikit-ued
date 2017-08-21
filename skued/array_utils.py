"""
Array utility functions
"""

from itertools import repeat
from numpy import concatenate, swapaxes, hypot, arctan2, sin, cos

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
        reverse = [slice(None, None, -1)] * arr.ndim
    else:
        reverse = [slice(None, None, None)] * arr.ndim

        if isinstance(axes, int):
            axes = (axes,)
            
        for axis in axes:
            reverse[axis] = slice(None, None, -1)
    
    return arr[reverse]

def cart2polar(x, y):
    """ 
    Transform cartesian coordinates to polar coordinates

    Parameters
    ----------
    x, y : `~numpy.ndarray`
        Cartesian coordinates
    
    Returns
    -------
    r, t : `~numpy.ndarray`
        Radius and polar angle coordinates
    """
    return hypot(x,y), arctan2(y, x)

def polar2cart(r, t):
    """
    Transform cartesian coordinates to polar coordinates

    Parameters
    ----------
    r, t : `~numpy.ndarray`
        Radius and polar angle coordinates

    Returns
    -------
    x, y : `~numpy.ndarray`
        Cartesian coordinates
    """
    return r * cos(t), r * sin(t)