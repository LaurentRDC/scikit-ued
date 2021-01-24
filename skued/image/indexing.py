# -*- coding: utf-8 -*-
"""
Indexing single crystals 
------------------------
"""
from warnings import catch_warnings, simplefilter

import numpy as np
from scipy.ndimage import gaussian_filter, laplace, shift
from skimage.registration._masked_phase_cross_correlation import cross_correlate_masked
import skimage.filters as filters
from skimage.measure import label, regionprops
from skimage.morphology import binary_erosion, disk

from ..fft import with_skued_fft
from .center import autocenter


@with_skued_fft
def bragg_peaks(im, mask=None, center=None):
    """
    Extract the position of Bragg peaks in a single-crystal diffraction pattern.

    .. versionadded:: 2.1.3

    Parameters
    ----------
    im : ndarray, shape (N,M)
        Single-crystal diffraction pattern.
    mask : ndarray, shape (N,M), dtype bool, optional
        Mask that evaluates to `True` on valid pixels of ``im``.
    center : 2-tuple, optional
        Center of the diffraction pattern, in ``(row, col)`` format.
        If ``None``, the center will be determined via :func:`autocenter`.

    Returns
    -------
    peaks : list of 2-tuples
        List of coordinates ``[row, col]`` for every detected peak, sorted
        in order of how close they are to the center of the image.

    References
    ----------
    Liu, Lai Chung. Chemistry in Action: Making Molecular Movies with Ultrafast
    Electron Diffraction and Data Science, Chapter 2. Springer Nature, 2020.
    """
    if center is None:
        center = autocenter(im=im, mask=mask)

    im = np.array(im, copy=True, dtype=np.float)
    im -= im.min()

    with catch_warnings():
        simplefilter("ignore", category=RuntimeWarning)
        im /= gaussian_filter(input=im, sigma=min(im.shape) / 20, truncate=2)
    im = np.nan_to_num(im, copy=False)

    autocorr = cross_correlate_masked(arr1=im, arr2=im, m1=mask, m2=mask, mode="same")
    autocorr = np.abs(autocorr)

    # The regions of interest are defined on the labels made
    # from the autocorrelation of the image. The center of the autocorr
    # is the center of the array; we need to correct the offset.
    # This also allows to use the mask on properties derived
    # from the autocorr
    autocorr = shift(
        autocorr,
        shift=np.asarray(center) - np.array(im.shape) / 2,
        order=1,
        mode="nearest",
    )

    laplacian = -1 * laplace(autocorr)
    regions = (laplacian > 0.02) * mask

    # To prevent noise from looking like actual peaks,
    # we erode labels using a small selection area
    regions = binary_erosion(regions, selem=disk(2))

    labels = label(regions, return_num=False)
    props = regionprops(label_image=labels, intensity_image=im)
    # TODO: prune regions to remove regions that are very close to each other
    return sorted(
        [
            prop.weighted_centroid
            for prop in props
            if not np.any(np.isnan(prop.weighted_centroid))
        ],
        key=lambda p: np.linalg.norm(p - center),
    )
