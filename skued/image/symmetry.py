# -*- coding: utf-8 -*-
"""
Image manipulation involving symmetry
=====================================
"""
import numpy as np
import npstreams as ns
from skimage.transform import rotate

from ..array_utils import mirror


def nfold(im, mod, center=None, mask=None, fill_value=0.0):
    """
    Returns an images averaged according to n-fold rotational symmetry. This can be used to
    boost the signal-to-noise ratio on an image with known symmetry, e.g. a diffraction pattern.

    Parameters
    ----------
    im : array_like, ndim 2
        Image to be azimuthally-symmetrized.
    center : array_like, shape (2,) or None, optional
        Coordinates of the center (in pixels) in the format ``center=[col, row]``. If ``center=None``,
        the image is rotated around the center of the array, i.e. ``center=(cols / 2 - 0.5, rows / 2 - 0.5)``.
    mod : int
        Fold symmetry number. Valid numbers must be a divisor of 360.
    mask : `~numpy.ndarray` or None, optional
        Mask of `image`. The mask should evaluate to `True` (or 1) on valid pixels.
        If None (default), no mask is used.
    fill_value : float, optional
        In the case of a mask that overlaps with itself when rotationally averaged,
        the overlapping regions will be filled with this value.

    Returns
    -------
    out : `~numpy.ndarray`, dtype float
        Symmetrized image.

    Raises
    ------
    ValueError : If `mod` is not a divisor of 360 deg.
    """
    if 360 % mod:
        raise ValueError(
            f"{mod}-fold rotational symmetry is not valid (not a divisor of 360)."
        )
    angles = range(0, 360, int(360 / mod))

    im = np.array(im, copy=True)

    kwargs = {"center": center, "mode": "constant", "cval": 0, "preserve_range": True}

    if mask is None:
        return ns.average(rotate(im, angle, **kwargs) for angle in angles)

    # Use weights because edges of the pictures, which might be cropped by the rotation
    # should not count in the average
    wt = np.ones_like(mask, dtype=im.dtype)
    wt[np.logical_not(mask)] = 0

    weights = (rotate(wt, angle, **kwargs) for angle in angles)
    imgs = (rotate(im, angle, **kwargs) for angle in angles)

    avg = ns.average(imgs, weights=weights)

    # Mask may overlap with itself during symmetrization. At those points, the average
    # will be zero (because the weights are 0 there)
    # However, because users may want to change that value to `fill_value != 0`, we need
    # to know where is the overlap
    if fill_value != 0.0:
        invalid_pixels = np.logical_not(mask)
        masks = (rotate(invalid_pixels, angle, **kwargs) for angle in angles)
        overlap = ns.prod(masks).astype(bool)  # equivalent to logical_and
        avg[overlap] = fill_value

    return ns.nan_to_num(avg, fill_value=fill_value)


def reflection(im, angle, center=None, mask=None, fill_value=0.0):
    """
    Symmetrize an image according to a reflection plane.

    Parameters
    ----------
    im : array_like, ndim 2
        Image to be symmetrized.
    angle : float
        Angle (in degrees) of the line that defines the reflection plane. This angle
        increases counter-clockwise from the positive x-axis. Angles
        larger that 360 are mapped back to [0, 360). Note that ``angle`` and ``angle + 180``
        are equivalent.
    center : array_like, shape (2,) or None, optional
        Coordinates of the center (in pixels). in the format ``center=[col, row]``. If ``center=None``,
        the image is rotated around the center of the array, i.e. ``center=(cols / 2 - 0.5, rows / 2 - 0.5)``.
    mask : `~numpy.ndarray` or None, optional
        Mask of `image`. The mask should evaluate to `True` (or 1) on valid pixels.
        If None (default), no mask is used.
    fill_value : float, optional
        In the case of a mask that overlaps with itself when rotationally averaged,
        the overlapping regions will be filled with this value.

    Returns
    -------
    out : `~numpy.ndarray`, dtype float
        Symmetrized image.
    """
    angle = float(angle) % 360

    # Data-type must be float because of use of NaN
    im = np.array(im, dtype=float, copy=True)
    reflected = np.array(im, copy=True)  # reflected image

    if mask is None:
        mask = np.ones_like(im, dtype=bool)
    invalid_pixels = np.logical_not(mask)

    kwargs = {"center": center, "mode": "constant", "cval": 0, "preserve_range": True}

    # Rotate the 'reflected' image so that the reflection line is the x-axis
    # Flip the image along the y-axis
    # Rotate back to original orientation
    # FIXME: this will not work properly for images that are offcenter
    def refl(arr):
        arr = rotate(arr, -angle, **kwargs)
        arr = mirror(arr, axes=0)
        arr = rotate(arr, angle, **kwargs)
        return arr

    reflected = refl(reflected)
    invalid_pixels_r = refl(invalid_pixels).astype(bool)

    result = ns.average([im, reflected])
    result[invalid_pixels_r] = fill_value

    return result
