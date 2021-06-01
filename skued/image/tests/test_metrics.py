# -*- coding: utf-8 -*-


import numpy as np
from npstreams import last

from skued import (
    combine_masks,
    isnr,
    mask_from_collection,
    mask_image,
    snr_from_collection,
    triml,
    trimr,
)


def test_snr_from_collection_trivial():
    """Test snr_from_collection on series of identical images"""
    images = [np.ones((64, 64)) for _ in range(10)]
    snr = snr_from_collection(images)

    assert np.allclose(snr, np.zeros_like(snr))


def test_snr_from_collection_correctess():
    """Test that snr_from_collection is equivalent to np.mean / np.std"""
    images = [np.random.random((64, 64)) for _ in range(10)]
    stack = np.dstack(images)

    from_numpy = np.mean(stack, axis=2) / np.std(stack, axis=2)
    from_skued = snr_from_collection(images)

    assert np.allclose(from_numpy, from_skued)


def test_isnr_trivial():
    """Test snr_from_collection on series of identical images"""
    images = [np.ones((64, 64)) for _ in range(10)]
    snr = last(isnr(images))

    assert np.allclose(snr, np.zeros_like(snr))


def test_isnr_correctess():
    """Test that snr_from_collection is equivalent to np.mean / np.std"""
    images = [np.random.random((64, 64)) for _ in range(10)]
    stack = np.dstack(images)

    from_numpy = np.mean(stack, axis=2) / np.std(stack, axis=2)
    from_skued = last(isnr(images))

    assert np.allclose(from_numpy, from_skued)


def test_mask_from_collection_trivial():
    """Test on set of images with value zero"""
    images = [np.zeros((64, 64)) for _ in range(5)]
    mask = mask_from_collection(images)

    assert images[0].shape == mask.shape
    assert not np.any(mask)


def test_mask_from_collection_intensity_threshold_no_lower_bound():
    """Test that intensity threshold is respected"""
    images = [np.ones((64, 64)) for _ in range(5)]
    images[2][32, 4] = 10
    mask = mask_from_collection(images, px_thresh=9)

    assert np.sum(mask) == 1  # only one pixels is masked
    assert mask[32, 4] == True


def test_mask_from_collection_intensity_threshold_with_lower_bound():
    """Test that intensity threshold is respected"""
    images = [np.ones((64, 64)) for _ in range(5)]
    images[2][32, 4] = -10
    mask = mask_from_collection(images, px_thresh=(0, np.inf))

    assert np.sum(mask) == 1  # only one pixels is masked
    assert mask[32, 4] == True


def test_mask_from_collection_std_threshold():
    """Test that std threshold is respected"""
    images = [np.ones((64, 64)) for _ in range(5)]
    images[0][5, 12] = 1000
    mask = mask_from_collection(images, px_thresh=10000, std_thresh=1)

    assert np.sum(mask) == 1  # only one pixels is masked
    assert mask[5, 12] == True


def test_mask_from_collection_single_image():
    """Test that mask_from_collection works even if input
    is a single array"""
    images = np.zeros((64, 64))
    mask = mask_from_collection(images)

    assert not np.any(mask)


def test_combine_masks_trivial():
    masks = tuple([np.zeros((64, 64)) for _ in range(5)])
    combined = combine_masks(*masks)

    assert not np.any(combined)


def test_combine_masks_single_element():
    masks = tuple([np.ones((64, 64)) for _ in range(5)])
    masks[0][4, 6] = False
    combined = combine_masks(*masks)

    assert np.sum(np.logical_not(combined)) == 1
    assert not combined[4, 6]


def test_mask_image_trivial():
    mask = np.ones((64, 64), dtype=bool)
    image = np.random.random((64, 64))
    masked = mask_image(image, mask)

    assert np.allclose(image, masked)


def test_mask_image_no_copy():
    """Test that mask_image can work in-place"""
    mask = np.random.randint(0, 1, size=(64, 64), dtype=bool)
    image = np.random.random((64, 64))
    masked = mask_image(image, mask, copy=False)

    assert image is masked


def test_trim_left_trivial():
    """Test that nothing is trimmed when percentile is zero"""
    array = np.arange(0, 101)
    trimmed = triml(array, percentile=0)
    assert np.allclose(array, trimmed)


def test_trim_left_all_axes():
    """Test that trimming is working"""
    array = np.arange(0, 101)  # [0, 1, 2, ..., 100]
    trimmed = triml(array, percentile=20, fill_value=100)
    assert np.all(trimmed >= 20)


def test_trim_left_fill_value():
    """Test that fill_value is indeed introduced"""
    array = np.arange(0, 101, dtype=float)  # [0, 1, 2, ..., 100]
    trimmed = triml(array, percentile=20, fill_value=np.nan)
    assert np.any(np.isnan(trimmed))


def test_trim_right_trivial():
    """Test that nothing is trimmed when percentile is 100"""
    array = np.arange(0, 101)
    trimmed = trimr(array, percentile=100)
    assert np.allclose(array, trimmed)


def test_trim_right_trimming():
    """Test that trimming is working"""
    array = np.arange(0, 101)  # [0, 1, 2, ..., 100]
    trimmed = trimr(array, percentile=20, fill_value=0)
    assert np.all(trimmed <= 20)


def test_trim_right_fill_value():
    """Test that fill_value is indeed introduced"""
    array = np.arange(0, 101, dtype=float)  # [0, 1, 2, ..., 100]
    trimmed = trimr(array, percentile=20, fill_value=np.nan)
    assert np.any(np.isnan(trimmed))
