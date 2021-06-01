# -*- coding: utf-8 -*-


import numpy as np

from skued.baseline.algorithms import (
    _dt_approx_rec,
    _dwt_approx_rec,
    _dwt_approx_rec2,
    baseline_dt,
    baseline_dwt,
)

np.random.seed(23)


def test_zero_level():
    """baseline computed at the zero-th level should not affect the array"""
    array = np.random.random(size=(101,))
    assert np.allclose(array, baseline_dwt(array, max_iter=5, level=0))


def test_approx_rec():
    """Test that the underlying _dwt_approx_rec function is working properly"""

    arr = np.random.random(size=(102,))
    rec_arr = _dwt_approx_rec(arr, level=2, wavelet="db6", mode="constant", axis=-1)
    assert rec_arr.shape == arr.shape

    arr2 = np.random.random(size=(101, 104))
    rec_arr2 = _dwt_approx_rec(arr2, level=2, wavelet="sym4", mode="constant", axis=1)
    assert rec_arr2.shape == arr2.shape

    arr2 = np.random.random(size=(102, 104))
    rec_arr2 = _dwt_approx_rec2(
        arr2, level=2, wavelet="sym6", mode="constant", axis=(-2, -1)
    )
    assert rec_arr2.shape == arr2.shape

    arr2 = np.random.random(size=(102, 94, 25))
    rec_arr2 = _dwt_approx_rec2(
        arr2, level=2, wavelet="sym6", mode="constant", axis=(0, 1)
    )
    assert rec_arr2.shape == arr2.shape


def test_trivial_case_1d():
    """The baseline for a 1d array of zeros should be zeros"""
    arr = np.zeros(shape=(101,))
    assert np.allclose(arr, baseline_dwt(arr, max_iter=10, level=None))


def test_trivial_case_2d():
    """The baseline for a 2d array of zeros should be zeros"""
    arr = np.zeros((21, 31, 5))
    assert np.allclose(arr, baseline_dwt(arr, max_iter=10, level=None, axis=(0, 1)))


def test_1d_along_axis():
    """Test that iterating over array rows and baseline_dwt along axis are equivalent"""
    block = np.random.random(size=(21, 51))

    # Iterate over rows
    baseline_iterated = np.empty_like(block)
    for index, row in enumerate(block):
        baseline_iterated[index, :] = baseline_dwt(row, max_iter=50)

        # along axis
    baseline_axis = baseline_dwt(block, max_iter=50, axis=1)

    assert np.allclose(baseline_axis, baseline_iterated)


def test_zero_level():
    """baseline computed at the zero-th level should not affect the array"""
    array = np.random.random(size=(101,))
    assert np.allclose(array, baseline_dt(array, max_iter=5, level=0))


def test_approx_rec():
    """Test that the underlying _dwt_approx_rec function is working properly"""

    arr = np.random.random(size=(102,))
    rec_arr = _dt_approx_rec(
        arr,
        level=2,
        first_stage="db1",
        wavelet="qshift3",
        mode="smooth",
        axis=-1,
    )
    assert rec_arr.shape == arr.shape

    arr2 = np.random.random(size=(21, 52))
    rec_arr2 = _dt_approx_rec(
        arr2,
        level=2,
        first_stage="db1",
        wavelet="qshift3",
        mode="smooth",
        axis=1,
    )
    assert rec_arr2.shape == arr2.shape


def test_trivial_case():
    """The baseline for an array of zeros should be zeros"""
    arr = np.zeros(shape=(101,))
    assert np.allclose(arr, baseline_dt(arr, max_iter=10, level=None))


def test_positive_baseline():
    """Test that the baseline is never negative"""
    arr = 10 * np.random.random(size=(128,))
    baseline = baseline_dt(arr, max_iter=10)

    assert np.all(np.greater_equal(baseline, 0))


def test_baseline_limit():
    """Test that the baseline is never more than the original signal"""
    arr = 10 * np.random.random(size=(128,))
    baseline = baseline_dt(arr, max_iter=10)

    assert np.all(np.greater_equal(arr, baseline))


def test_final_shape():
    """Test that baseline_dt returns an array of the same shape as input"""

    arr = np.random.random(size=(101,))
    b = baseline_dt(arr, max_iter=10, level=None)
    assert arr.shape == b.shape

    arr = np.random.random(size=(100,))
    b = baseline_dt(arr, max_iter=10, level=None)
    assert arr.shape == b.shape

    arr = np.random.random(size=(101, 113))
    b = baseline_dt(arr, max_iter=10, level=None, axis=0)
    assert arr.shape == b.shape

    arr = np.random.random(size=(102, 104))
    b = baseline_dt(arr, max_iter=10, level=None, axis=1)
    assert arr.shape == b.shape


def test_2d_along_axis():
    """Test that iterating over array rows and baseline_dt along axis are equivalent"""
    block = np.random.random(size=(21, 51))

    # Iterate over rows
    baseline_iterated = np.empty_like(block)
    for index, row in enumerate(block):
        baseline_iterated[index, :] = baseline_dt(row, max_iter=50)

        # along axis
    baseline_axis = baseline_dt(block, max_iter=50, axis=1)

    assert np.allclose(baseline_axis, baseline_iterated)


def test_background_regions():
    """Test that background_regions is used correctly along certain axes"""
    block = np.random.random(size=(21, 51))

    baseline_params = {"max_iter": 50, "background_regions": [(..., range(0, 3))]}

    # Iterate over rows
    baseline_iterated = np.empty_like(block)
    for index, row in enumerate(block):
        baseline_iterated[index, :] = baseline_dt(row, **baseline_params)

        # along axis
    baseline_axis = baseline_dt(block, axis=1, **baseline_params)

    assert np.allclose(baseline_axis, baseline_iterated)
