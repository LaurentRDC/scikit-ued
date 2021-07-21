# -*- coding: utf-8 -*-


import numpy as np
import pywt
import pytest

from skued.baseline import (
    available_dt_filters,
    available_first_stage_filters,
    dt_first_stage,
    dt_max_level,
    dtcwt,
    idtcwt,
)
import pytest

np.random.seed(23)


@pytest.mark.parametrize("first_stage_wavelet", available_first_stage_filters())
def test_first_stage(first_stage_wavelet):
    """Test of perfect reconstruction of first stage wavelets."""
    array = np.sin(np.arange(0, 10, step=0.01))
    # Using waverec and wavedec instead of dwt and idwt because parameters
    # don't need as much parsing.
    assert np.allclose(
        array,
        pywt.waverec(pywt.wavedec(array, first_stage_wavelet), first_stage_wavelet),
    )


@pytest.mark.parametrize("first_stage_wavelet", available_first_stage_filters())
def test_first_stage_issue_36(first_stage_wavelet):
    """Test that first-stage wavelets are properly shifted. See Issue 36"""
    w1, w2 = dt_first_stage(first_stage_wavelet)
    # Reconstruction should be shifted back
    for f1, f2 in zip(w1.filter_bank[:2], w2.filter_bank[:2]):
        assert np.allclose(f1[1:-1], f2[2::])

    # Deconstruction should be shifted forward
    for f1, f2 in zip(w1.filter_bank[2::], w2.filter_bank[2::]):
        assert np.allclose(f1[1:-1], f2[0:-2])


def gen_input(n_dimensions):
    """Generate random array with the appropriate dimensions"""
    shape = {1: (100,), 2: (50, 50), 3: (10, 10, 10)}[n_dimensions]
    return np.random.random(size=shape)


multidim = pytest.mark.parametrize("n_dimensions", (1, 2, 3))


@multidim
def test_perfect_reconstruction_level_0(n_dimensions):
    """Test perfect reconstruction for a 0-level decomposition"""
    array = gen_input(n_dimensions)
    coeffs = dtcwt(data=array, level=0, first_stage="sym6", wavelet="qshift1")
    reconstructed = idtcwt(coeffs=coeffs, first_stage="sym6", wavelet="qshift1")
    assert np.allclose(array, reconstructed)


@multidim
def test_perfect_reconstruction_level_1(n_dimensions):
    """Test perfect reconstruction for a single decomposition level"""
    array = gen_input(n_dimensions)
    for first_stage in available_first_stage_filters():
        coeffs = dtcwt(data=array, level=1, first_stage=first_stage, wavelet="qshift1")
        reconstructed = idtcwt(
            coeffs=coeffs, first_stage=first_stage, wavelet="qshift1"
        )
        assert np.allclose(array, reconstructed)


@multidim
def test_perfect_reconstruction_multilevel(n_dimensions):
    """Test perfect reconstruction for all levels, for all first_stage wavelets, for all DT wavelets"""
    array = gen_input(n_dimensions)

    for first_stage in available_first_stage_filters():
        for wavelet in available_dt_filters():
            for level in range(
                1,
                dt_max_level(data=array, first_stage=first_stage, wavelet=wavelet),
            ):
                coeffs = dtcwt(
                    data=array,
                    level=level,
                    first_stage=first_stage,
                    wavelet=wavelet,
                )
                reconstructed = idtcwt(
                    coeffs=coeffs, first_stage=first_stage, wavelet=wavelet
                )
                assert np.allclose(array, reconstructed)


@multidim
def test_axis(n_dimensions):
    """Test perfect reconstruction along all axes"""
    array = gen_input(n_dimensions)
    for axis in range(0, array.ndim):
        coeffs = dtcwt(
            data=array,
            level=2,
            axis=axis,
            first_stage="sym6",
            wavelet="qshift1",
        )
        reconstructed = idtcwt(
            coeffs=coeffs, axis=axis, first_stage="sym6", wavelet="qshift1"
        )
        assert np.allclose(array, reconstructed)


@multidim
def test_axis_limits(n_dimensions):
    """Test that an exception is raised for an invalid 'axis' parameter"""
    array = gen_input(n_dimensions)
    with pytest.raises(ValueError):
        coeffs = dtcwt(
            data=array,
            level=1,
            axis=array.ndim,
            first_stage="sym6",
            wavelet="qshift1",
        )


@multidim
def test_even_length_along_axis(n_dimensions):
    """Test that an exception is raised when array is not even along transform axis"""
    with pytest.raises(ValueError):
        dtcwt(
            data=np.zeros((17, 8)),
            level=1,
            first_stage="sym6",
            wavelet="qshift1",
            axis=0,
        )
