# -*- coding: utf-8 -*-


import numpy as np

from skued import mad


def test_mad_trivial():
    """Test that the median absolute dev of an array of zeroes is zero"""
    arr = np.zeros((16,))
    assert np.allclose(mad(arr), 0)


def test_mad_integers():
    """Test that mad of an integer array is working as intended"""
    arr_int = np.random.randint(-15, 15, size=(8, 8))
    arr_flo = np.array(arr_int, copy=True, dtype=float)

    mad_int = mad(arr_int)
    mad_flo = mad(arr_flo)

    assert np.allclose(mad_int, mad_flo)


def test_mad_correctness():
    """Test that the mad() function works as expected"""
    arr = np.random.random(size=(32,))
    from_mad = mad(arr)
    explicit = np.abs(arr - np.median(arr))

    assert np.allclose(from_mad, explicit)


def test_mad_side_effects():
    """Test that input array is not modified by mad()"""
    arr = np.random.random(size=(32,))
    arr.setflags(write=False)
    mad(arr)
