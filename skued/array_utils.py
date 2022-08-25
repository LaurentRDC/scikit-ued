"""
Array utility functions
"""

from itertools import repeat
from warnings import warn

import numpy as np
from numpy.linalg import norm


def repeated_array(arr, num, axes=-1):
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
    ValueError : If num and axes are tuples of different lengths.
    """
    if not num:
        return arr

    if isinstance(num, int):
        num = (num,)
    if isinstance(axes, int):
        axes = (axes,)

    if len(num) != len(axes):
        raise ValueError("num and axes must have the same length")

    composite = np.concatenate(tuple(repeat(arr, times=num[0])), axis=axes[0])

    if len(num) > 1:
        for n, ax in zip(num[1:], axes[1:]):
            composite = np.concatenate(tuple(repeat(composite, times=n)), axis=ax)

    return composite


def complex_array(real, imag):
    """
    Combine two real ndarrays into a complex array.

    Parameters
    ----------
    real, imag : array_like
        Real and imaginary parts of a complex array.

    Returns
    -------
    complex : `~numpy.ndarray`
        Complex array.
    """
    real, imag = np.asfarray(real), np.asfarray(imag)
    comp = real.astype(complex)
    comp += 1j * imag
    return comp


def mirror(arr, axes=None):
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
    out : `~numpy.ndarray`
        Mirrored array.
    """
    if axes is None:
        reverse = [slice(None, None, -1)] * arr.ndim
    else:
        reverse = [slice(None, None, None)] * arr.ndim

        if isinstance(axes, int):
            axes = (axes,)

        for axis in axes:
            reverse[axis] = slice(None, None, -1)

    return arr[tuple(reverse)]


def cart2polar(x, y):
    """
    Transform cartesian coordinates to polar coordinates.

    Parameters
    ----------
    x, y : `~numpy.ndarray`
        Cartesian coordinates.

    Returns
    -------
    r, t : `~numpy.ndarray`
        Radius and polar angle coordinates. Polar angle coordinates
        are in radians.
    """
    return np.hypot(x, y), np.arctan2(y, x)


def polar2cart(r, t):
    """
    Transform polar coordinates to cartesian coordinates.

    Parameters
    ----------
    r, t : `~numpy.ndarray`
        Radius and polar angle coordinates. Angles
        are assumed to be radians.

    Returns
    -------
    x, y : `~numpy.ndarray`
        Cartesian coordinates
    """
    return r * np.cos(t), r * np.sin(t)


def spherical2cart(r, p, t):
    """
    Transform spherical coordinates into cartesian coordinates.

    Parameters
    ----------
    r : `~numpy.ndarray`
        Radial coordinate.
    p : `~numpy.ndarray`
        Polar angle coordinate in radians.
    t : `~numpy.ndarray`
        Azimuthal coordinate in radians.

    Returns
    -------
    x, y, z : `~numpy.ndarray`
        Cartesian coordinates.
    """
    x = r * np.sin(t) * np.cos(p)
    y = r * np.sin(t) * np.sin(p)
    z = r * np.cos(t)
    return x, y, z


def cart2spherical(x, y, z):
    """
    Transform cartesian coordinates into spherical coordinates .

    Parameters
    ----------
    x, y, z : `~numpy.ndarray`
        Cartesian coordinates

    Returns
    -------
    r : `~numpy.ndarray`
        Radial coordinate.
    p : `~numpy.ndarray`
        Polar angle coordinate in radians.
    t : `~numpy.ndarray`
        Azimuthal coordinate in radians.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    p = np.arctan2(y, x)
    t = np.arccos(z / r)
    return r, p, t


def plane_mesh(v1, v2, x1, x2=None, origin=(0, 0, 0)):
    """
    Generate a spatial mesh for a plane defined by two vectors.

    Parameters
    ----------
    v1, v2 : `~numpy.ndarray`, shape (3,)
        Basis vector of the plane. A warning is raised if
        ``v1`` and ``v2`` are not orthogonal.
    x1, x2 : iterable, shape (N,)
        1-D arrays representing the coordinates on the grid, along basis
        vectors ``v1`` and ``v2`` respectively. If ``x2`` is not specified,
        ``x1`` and ``x2`` are taken to be the same
    origin : iterable, shape (3,), optional
        Plane mesh will be generated with respect to this origin.

    Returns
    -------
    x, y, z : `~numpy.ndarray`, ndim 2
        Mesh arrays for the coordinate of the plane.
    """
    v1, v2 = v1 / norm(v1), v2 / norm(v2)

    if x2 is None:
        x2 = np.array(x1)

    if np.dot(v1, v2) != 0:
        warn("Plane basis vectors are not orthogonal", RuntimeWarning)

    along_v1, along_v2 = np.meshgrid(x1, x2, indexing="ij")
    xx, yy, zz = tuple(along_v1 * v1[i] + along_v2 * v2[i] for i in range(3))

    ox, oy, oz = tuple(origin)
    return (xx + ox, yy + oy, zz + oz)
