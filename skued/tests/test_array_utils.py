# -*- coding: utf-8 -*-


import numpy as np

from skued import (
    cart2polar,
    cart2spherical,
    complex_array,
    mirror,
    plane_mesh,
    polar2cart,
    repeated_array,
    spherical2cart,
)

np.random.seed(23)


def test_repeated_array_trivial():
    """Test repeated_array of 0 copies"""
    arr = np.random.random(size=(4, 5))
    composite = repeated_array(arr, num=0, axes=0)
    assert np.allclose(composite, arr)


def test_repeated_array_single_axis():
    """Test repeating the array over a single axis"""
    arr = np.random.random(size=(4, 5))
    composite = repeated_array(arr, num=3, axes=1)
    expected_new_shape = (arr.shape[0], arr.shape[1] * 3)
    assert composite.shape == expected_new_shape


def test_repeated_array_multiple_axes():
    """Test repeating an array over multiple axes in all possible orders"""
    arr = np.random.random(size=(4, 5))
    composite = repeated_array(arr, num=(3, 2), axes=(0, 1))
    expected_new_shape = (arr.shape[0] * 3, arr.shape[1] * 2)
    assert composite.shape == expected_new_shape

    composite = repeated_array(arr, num=(2, 3), axes=(1, 0))
    expected_new_shape = (arr.shape[0] * 3, arr.shape[1] * 2)
    assert composite.shape == expected_new_shape


def test_complex_array_floats():
    """Test that two floating arrays are cast correctly"""
    real, imag = np.empty((3, 4), dtype=float), np.empty((3, 4), dtype=float)
    assert complex_array(real, imag).dtype == complex


def test_complex_array_non_floats():
    """Test that two integer arrays are cast correctly"""
    real, imag = np.empty((3, 4), dtype=np.int16), np.empty((3, 4), dtype=np.int8)
    assert complex_array(real, imag).dtype == complex


def test_complex_array_results():
    """test that ``complex_array`` returns appropriate results"""
    arr1 = np.random.random((4, 5))
    arr2 = np.random.random((4, 5))

    from_complex_array = complex_array(arr1, arr2)
    by_hand = arr1 + 1j * arr2
    assert from_complex_array.dtype == by_hand.dtype
    assert np.allclose(from_complex_array, by_hand)


def test_mirror_1D():
    """Test mirror() on a 1D array"""
    arr = np.zeros((16,), dtype=float)
    arr[15] = 1
    assert np.allclose(arr[::-1], mirror(arr))


def test_mirror_2D_all_axes():
    """Test mirror() on a 2D array for all axes"""
    arr = np.zeros((16, 16), dtype=float)
    arr[15, 3] = 1
    assert np.allclose(arr[::-1, ::-1], mirror(arr))


def test_mirror_2D_one_axis():
    """Test mirror() on a 2D array for one axis"""
    arr = np.zeros((16, 16), dtype=float)
    arr[15, 3] = 1
    assert np.allclose(arr[:, ::-1], mirror(arr, axes=1))
    assert np.allclose(arr[::-1, :], mirror(arr, axes=0))


def test_cart2polar_back_and_forth():
    """Test that cart2polar and polar2cart are reciprocal"""
    x = np.random.random(size=(16, 8))
    y = np.random.random(size=(16, 8))

    r, t = cart2polar(x, y)

    xp, yp = polar2cart(r, t)
    assert np.allclose(x, xp)
    assert np.allclose(y, yp)


def test_spherical2cart_back_and_forth():
    """Test that cart2polar and polar2cart are reciprocal"""
    x = np.random.random(size=(16, 8))
    y = np.random.random(size=(16, 8))
    z = np.random.random(size=(16, 8))

    r, p, t = cart2spherical(x, y, z)

    xp, yp, zp = spherical2cart(r, p, t)
    assert np.allclose(x, xp)
    assert np.allclose(y, yp)
    assert np.allclose(z, zp)


def test_plane_mesh_shape():
    """Test that shape is as expected"""
    extent1 = np.linspace(0, 10, num=64)
    extent2 = np.linspace(0, 10, num=128)
    v1, v2, _ = np.eye(3)

    for arr in plane_mesh(v1, v2, extent1, extent2):
        assert arr.shape == (64, 128)


def test_plane_mesh_origin():
    """Test that plane_mesh is generated from origin"""
    extent1 = np.linspace(0, 10, num=64)
    extent2 = np.linspace(0, 10, num=128)
    v1, v2, _ = np.eye(3)

    for arr in plane_mesh(v1, v2, extent1, extent2, origin=(-4, -4, -4)):
        assert arr.min() == -4
