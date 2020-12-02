# -*- coding: utf-8 -*-
"""
Image manipulation involving symmetry
=====================================
"""
import numpy as np
from npstreams import average, nan_to_num
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
        Coordinates of the center (in pixels). If ``center=None``, the image is rotated around
        its center, i.e. ``center=(rows / 2 - 0.5, cols / 2 - 0.5)``.
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

    # Data-type must be float because of use of NaN
    im = np.array(im, dtype=np.float, copy=True)

    if mask is None:
        mask = np.ones_like(im, dtype=np.uint8)

    kwargs = {"center": center, "mode": "constant", "cval": 0, "preserve_range": True}

    # Use weights because edges of the pictures, which might be cropped by the rotation
    # should not count in the average
    wt = np.ones_like(mask, dtype=np.float)
    wt[np.logical_not(mask)] = np.nan

    weights = (rotate(wt, angle, **kwargs) for angle in angles)
    rotated = (rotate(im, angle, **kwargs) for angle in angles)

    avg = average(rotated, weights=weights, ignore_nan=True)
    return nan_to_num(avg, fill_value, copy=False)


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
        Coordinates of the center (in pixels). If ``center=None``, the image is rotated around
        its center, i.e. ``center=(rows / 2 - 0.5, cols / 2 - 0.5)``.
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
    im = np.array(im, dtype=np.float, copy=True)
    reflected = np.array(im, copy=True)  # reflected image

    if mask is None:
        mask = np.ones_like(im, dtype=np.bool)
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
    invalid_pixels_r = refl(invalid_pixels).astype(np.bool)

    result = average([im, reflected])
    result[invalid_pixels_r] = fill_value

    return result
