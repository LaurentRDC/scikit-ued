# -*- coding: utf-8 -*-

from random import randint

import numpy as np
from scipy import ndimage as ndi
from skimage import data
from skimage.draw import ellipse
from skimage.filters import gaussian
import skimage.data as data

from skued import align, ialign, itrack_peak


np.random.seed(23)


def circle_image(shape, center, radii, intensities):
    """Creates an image with circle or thickness 2"""
    im = np.zeros(shape=shape, dtype=float)
    xx, yy = np.ogrid[0 : shape[0], 0 : shape[1]]
    xx, yy = xx - center[0], yy - center[1]
    for radius, intensity in zip(radii, intensities):
        rr = np.sqrt(xx**2 + yy**2)
        im[np.logical_and(rr < radius + 1, rr > radius - 1)] = intensity

    im[:] = gaussian(im, 5)
    return im


def camera():
    return data.camera()[::2, ::2]  # decimate by 2 for faster tests


def test_ialign_trivial():
    """Test alignment of identical images"""
    aligned = tuple(ialign([camera() for _ in range(5)]))

    assert len(aligned) == 5
    assert camera().shape == aligned[0].shape


def test_ialign_misaligned_canned_images():
    """shift images from skimage.data by entire pixels.
    We don't expect perfect alignment."""
    reference = camera()
    shifts = [(randint(-4, 0), randint(0, 4)) for _ in range(5)]

    mask = np.zeros_like(reference, dtype=bool)
    rr1, cc1 = ellipse(129, 127, r_radius=63, c_radius=50, shape=reference.shape)
    mask[rr1, cc1] = True

    misaligned = (ndi.shift(camera(), shift=s) for s in shifts)

    for aligned, (sx, sy) in zip(
        ialign(misaligned, reference=reference, mask=mask), shifts
    ):
        assert np.allclose(reference[sx::, 0:-sy], aligned[sx::, 0:-sy])


def test_align_no_side_effects():
    """Test that aligned images are not modified in-place"""
    im = np.array(camera()[0:64, 0:64])
    im.setflags(write=False)
    aligned = align(im, reference=im, fill_value=np.nan)
    assert im.dtype == camera().dtype


def test_align_no_mask():
    """Test that alignment of images with no masks works"""
    reference = camera()
    image = ndi.shift(camera(), shift=(-7, 12))
    aligned = align(image, reference=reference)

    assert np.allclose(reference[7::, 0:-12], aligned[7::, 0:-12])


def test_align_with_mask():
    """Test that alignment of images with no masks works"""
    reference = camera()
    image = ndi.shift(camera(), shift=(-7, 12))

    mask = np.ones_like(reference, dtype=bool)
    mask[75:100, 50:100] = False

    aligned = align(image, reference=reference, mask=mask)

    assert np.allclose(reference[7::, 0:-12], aligned[7::, 0:-12])


def test_itrack_peak_trivial():
    """Test that shift is identically zero for images that are identical"""
    # Array prototype is just zeros
    # with a 'peak' in the center
    prototype = np.zeros(shape=(17, 17))
    prototype[9, 9] = 10
    images = [np.array(prototype) for _ in range(20)]
    shifts = itrack_peak(images, row_slice=np.s_[:], col_slice=np.s_[:])

    for shift in shifts:
        assert np.allclose(shift, (0.0, 0.0))


def test_itrack_peak_length():
    """Test that shifts yielded by itrack_peak are as numerous as the number of input pictures"""
    images = [np.random.random(size=(4, 4)) for _ in range(20)]
    shifts = list(itrack_peak(images, row_slice=np.s_[:], col_slice=np.s_[:]))

    assert len(shifts) == len(images)
