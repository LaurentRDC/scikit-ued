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
def bragg_peaks(im, mask=None, center=None, min_dist=None):
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
    min_dist : float or None, optional
        Minimum distance between Bragg peaks (in pixel coordinates). Peaks that are closer
        than this distance will be considered the same peak, and only one of them will
        be returned. If `None` (default), the minimum distance is guessed based on the
        image size.

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

    im = np.array(im, copy=True, dtype=float)
    im -= im.min()

    with catch_warnings():
        simplefilter("ignore", category=RuntimeWarning)
        im /= gaussian_filter(input=im, sigma=min(im.shape) / 20, truncate=2)
    im = np.nan_to_num(im, copy=False)

    autocorr = np.abs(
        cross_correlate_masked(arr1=im, arr2=im, m1=mask, m2=mask, mode="same")
    )

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
    threshold = filters.threshold_triangle(laplacian)
    regions = (laplacian > threshold) * mask

    # To prevent noise from looking like actual peaks,
    # we erode labels using a small selection area
    regions = binary_erosion(regions, selem=disk(2))

    labels = label(regions, return_num=False)
    props = regionprops(label_image=labels, intensity_image=im)
    candidates = [
        prop for prop in props if not np.any(np.isnan(prop.weighted_centroid))
    ]

    # Some regions are very close to each other; we prune them!
    if min_dist is None:
        min_dist = min(im.shape) / 100

    peaks = list()
    for prop in candidates:
        pos = np.asarray(prop.weighted_centroid)
        if any((np.linalg.norm(peak - pos) < min_dist) for peak in peaks):
            continue
        else:
            peaks.append(pos)

    return sorted(peaks, key=lambda p: np.linalg.norm(p - center))
