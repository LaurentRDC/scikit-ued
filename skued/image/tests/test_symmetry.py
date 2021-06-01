# -*- coding: utf-8 -*-

from warnings import catch_warnings, simplefilter

import numpy as np

from skued import nfold, reflection
import pytest

np.random.seed(23)


def test_nfold_trivial():
    """Test nfold_symmetry averaging on trivial array"""
    im = np.zeros((256, 256))
    rot = nfold(im, mod=3)
    assert np.allclose(rot, im)


def test_nfold_valid_mod():
    """Test the the N-fold symmetry argument is valid"""
    im = np.empty((128, 128))
    with pytest.raises(ValueError):
        nfold(im, mod=1.7)


def test_nfold_mask():
    """Test that nfold() works correctly with a mask"""
    im = np.zeros((128, 128), dtype=int)
    mask = np.ones_like(im, dtype=bool)

    with catch_warnings():
        simplefilter("ignore")
        im[0:20] = 1
        mask[0:20] = False

        rot = nfold(im, mod=2, mask=mask)
        assert np.allclose(rot, np.zeros_like(rot))


def test_nfold_mask_with_overlapping_symmetry():
    """Test that nfold() works correctly with a mask which overlaps with it when rotated."""
    im = np.zeros((128, 128), dtype=int)
    mask = np.ones_like(im, dtype=bool)

    with catch_warnings():
        simplefilter("ignore")
        im[64 - 4 : 64 + 4, :] = 1
        mask[64 - 5 : 64 + 5, :] = False  # mask slightly larger

        rot = nfold(im, mod=6, mask=mask)
        assert np.allclose(rot, np.zeros_like(rot))


def test_nfold_mask_with_overlapping_symmetry_fill_value():
    """Test that nfold() works correctly with a mask which overlaps with it when rotated, when `fill_value` is not 0."""
    im = np.zeros((128, 128), dtype=int)
    mask = np.ones_like(im, dtype=bool)

    with catch_warnings():
        simplefilter("ignore")
        im[64 - 4 : 64 + 4, :] = 1
        mask[64 - 5 : 64 + 5, :] = False  # mask slightly larger

        rot = nfold(im, mod=6, mask=mask, fill_value=500)
        assert rot[64, 64] == 500


def test_nfold_no_side_effects():
    """Test that nfold() does not modify the input image and mask"""
    im = np.empty((128, 128), dtype=float)
    mask = np.zeros_like(im, dtype=bool)

    im.setflags(write=False)
    mask.setflags(write=False)

    with catch_warnings():
        simplefilter("ignore")
        rot = nfold(im, center=(67, 93), mod=3, mask=mask)


def test_nfold_fill_value():
    """Test that the fill_value parameter of nfold() is working correctly"""
    im = 1000 * np.random.random(size=(256, 256))
    mask = np.random.choice([True, False], size=im.shape)

    with catch_warnings():
        simplefilter("ignore")
        rot = nfold(im, center=(100, 150), mod=5, mask=mask, fill_value=np.nan)

    assert np.any(np.isnan(rot))

    with catch_warnings():
        simplefilter("ignore")
        rot = nfold(im, center=(100, 150), mod=5, mask=mask, fill_value=0.0)

        assert not np.any(np.isnan(rot))


def test_nfold_output_range():
    """Test that nfold() does not modify the value range"""
    im = 1000 * np.random.random(size=(256, 256))
    mask = np.random.choice([True, False], size=im.shape)

    with catch_warnings():
        simplefilter("ignore")
        rot = nfold(im, center=(100, 150), mod=5, mask=mask)

    assert rot.max() <= im.max()
    # In the case of a mask that overlaps with it when rotated,
    # the average will be zero due to nan_to_num
    assert rot.min() >= min(im.min(), 0)


def test_nfold_mod_1():
    """Test that nfold(mod = 1) returns an unchanged image, except
    perhaps for a cast to float"""
    im = 1000 * np.random.random(size=(256, 256))
    rot = nfold(im, mod=1)
    assert np.allclose(im, rot)


def test_reflection_trivial():
    """Test the reflection symmetry of an image of zeroes"""
    im = np.zeros((256, 256))
    ref = reflection(im, angle=-15.3)
    assert np.allclose(ref, im)


def test_reflection_no_side_effects():
    """Test that reflection() does not modify the input image and mask"""
    im = np.empty((128, 128), dtype=float)
    mask = np.zeros_like(im, dtype=bool)

    im.setflags(write=False)
    mask.setflags(write=False)

    rot = reflection(im, center=(67, 93), angle=35, mask=mask)


def test_reflection_output_range():
    """Test that reflection() does not modify the value range"""
    im = 1000 * np.random.random(size=(256, 256))
    mask = np.random.choice([True, False], size=im.shape)

    with catch_warnings():
        simplefilter("ignore")
        rot = reflection(im, center=(100, 150), angle=5, mask=mask)

    assert rot.max() <= im.max()
    # In the case of a mask that overlaps with it when rotated,
    # the average will be zero due to nan_to_num
    assert rot.min() >= min(im.min(), 0)


def test_reflection_correctness_angle0():
    """Test that reflection() correctly symmetrizes around the x-axis."""

    im = np.zeros((256, 256), dtype=float)
    im[0:10, :] = 1

    reflected = reflection(im, angle=0)

    expected = np.array(im, copy=True)
    expected[0:10, :] = 0.5
    expected[246:256, :] = 0.5

    assert np.allclose(expected, reflected)


def test_reflection_angle_vs_angle_plus_180():
    """Test that the result of reflection() is the same for any angle and angle + 180"""
    im = np.random.random((256, 256))
    reflected1 = reflection(im, angle=15)
    reflected2 = reflection(im, angle=195)  # + 180

    assert np.allclose(reflected1, reflected2)


def test_reflection_correctness_angle90():
    """Test that reflection() correctly symmetrizes around the y-axis."""

    im = np.zeros((256, 256), dtype=float)
    im[:, 0:10] = 1

    reflected = reflection(im, angle=90)

    expected = np.array(im, copy=True)
    expected[:, 0:10] = 0.5
    expected[:, 246:256] = 0.5

    assert np.allclose(expected, reflected)
