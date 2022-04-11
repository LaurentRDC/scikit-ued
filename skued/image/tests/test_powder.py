# -*- coding: utf-8 -*-


import numpy as np
from skimage.filters import gaussian

from skued import azimuthal_average

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


def test_azimuthal_average_trivial_array():
    """Test azimuthal_average on an array of zeroes"""
    image = np.zeros(shape=(256, 256), dtype=float)
    center = (image.shape[0] / 2, image.shape[1] / 2)

    radius, intensity = azimuthal_average(image, center)

    assert intensity.sum() == 0
    assert intensity.shape == radius.shape


def test_azimuthal_average_ring():
    """Test azimuthal_average on an image with a wide ring"""
    image = np.zeros(shape=(256, 256), dtype=float)
    center = (image.shape[0] / 2, image.shape[1] / 2)
    xc, yc = center

    # Create an image with a wide ring
    extent = np.arange(0, image.shape[0])
    xx, yy = np.meshgrid(extent, extent)
    rr = np.sqrt((xx - xc) ** 2 + (yy - yc) ** 2)
    image[np.logical_and(24 < rr, rr < 26)] = 1

    radius, intensity = azimuthal_average(image, center)
    assert intensity.max() == image.max()
    assert radius.shape == intensity.shape


def test_azimuthal_average_angular_bounds():
    """Test azimuthal_average with a restrictive angular_bounds argument"""
    image = np.zeros(shape=(256, 256), dtype=float)
    center = (image.shape[0] / 2, image.shape[1] / 2)
    xc, yc = center

    # Create an image with a wide ring
    extent = np.arange(0, image.shape[0])
    xx, yy = np.meshgrid(extent, extent)
    rr = np.sqrt((xx - xc) ** 2 + (yy - yc) ** 2)
    angles = np.rad2deg(np.arctan2(yy - yc, xx - xc)) + 180
    image[np.logical_and(0 <= angles, angles <= 60)] = 1

    radius, intensity = azimuthal_average(image, center, angular_bounds=None)
    r360, int360 = azimuthal_average(image, center, angular_bounds=(0, 360))
    assert np.allclose(intensity, int360)

    radius, intensity = azimuthal_average(image, center, angular_bounds=(0, 60))
    assert np.allclose(intensity, np.ones_like(intensity))

    radius, intensity = azimuthal_average(image, center, angular_bounds=(15, 75))
    assert not np.all(intensity < np.ones_like(intensity))

    radius, intensity = azimuthal_average(image, center, angular_bounds=(60, 360))
    assert np.allclose(intensity, np.zeros_like(intensity))

    radius, intensity = azimuthal_average(
        image, center, angular_bounds=(60 + 360, 360 + 360)
    )
    assert np.allclose(intensity, np.zeros_like(intensity))


def test_azimuthal_average_ring_with_mask():
    """Test azimuthal_average on an image with a wide ring"""
    image = np.zeros(shape=(256, 256), dtype=float)
    center = (image.shape[0] / 2, image.shape[1] / 2)
    xc, yc = center

    mask = np.ones_like(image, dtype=bool)
    mask[120:140, 0:140] = False

    # Create an image with a wide ring
    extent = np.arange(0, image.shape[0])
    xx, yy = np.meshgrid(extent, extent)
    rr = np.sqrt((xx - xc) ** 2 + (yy - yc) ** 2)
    image[np.logical_and(24 < rr, rr < 26)] = 1
    # Add some ridiculously high value under the mask
    # to see if it will be taken intou account
    image[128, 128] = 10000

    radius, intensity = azimuthal_average(image, center, mask=mask, trim=False)

    # The maximum value outside of the mask area was set to 1
    assert intensity.max() == 1
    assert radius.shape == intensity.shape


def test_azimuthal_average_trim_and_mask():
    """Test that regions that only have masks contributions are not present
    in the angular average"""
    image = np.ones(shape=(256, 256), dtype=float)
    center = (image.shape[0] / 2, image.shape[1] / 2)
    xc, yc = center

    # Create an image with a wide ring
    extent = np.arange(0, image.shape[0])
    xx, yy = np.meshgrid(extent, extent)
    rr = np.hypot(xx - xc, yy - yc)

    mask = np.ones_like(image, dtype=bool)
    mask[rr < 20] = False
    # image[rr < 20] = 0

    radius, intensity = azimuthal_average(image, center, mask=mask, trim=False)
    assert radius.min() == 0

    radius_trimmed, intensity_trimmed = azimuthal_average(
        image, center, mask=mask, trim=True
    )
    assert radius_trimmed.min() == 20


def test_azimuthal_average_mask_and_nan():
    """Test that azimuthal_average with masks does not yield NaNs. This can happen for large masks."""
    image = np.ones(shape=(256, 256), dtype=np.int16)
    mask = np.zeros_like(image, dtype=bool)
    mask[100:156, 100:156] = True
    _, av = azimuthal_average(image, center=(128, 128), mask=mask, trim=False)

    assert not np.any(np.isnan(av))
