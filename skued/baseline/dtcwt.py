# -*- coding: utf-8 -*-
""" 
Dual-tree complex wavelet transform
===================================
"""
from collections import deque
from functools import lru_cache
from itertools import cycle
from pathlib import Path

import numpy as np
from pywt import Wavelet, dwt, dwt_max_level, idwt, wavelist

DATADIR = Path(__file__).parent / "data"


def available_dt_filters():
    """
    Iterable of available wavelet filters compatible with the dual-tree complex wavelet transform.

    Returns
    -------
    wavelets : iterable
        List of sorted string names. The wavelet numerical values can be
        retrieved from the :func:`dualtree_wavelet` function.
    """
    return sorted(["qshift1", "qshift2", "qshift3", "qshift4", "qshift5", "qshift6"])


def available_first_stage_filters():
    """
    Iterable of available wavelet filters compatible with the
    first stage of dual-tree complex wavelent transform.

    Returns
    -------
    wavelets : iterable
        List of sorted string names. The wavelet numerical values can be
        retrieved from the :func:`dt_first_stage` function.
    """
    return sorted(filter(lambda name: name != "dmey", wavelist(kind="discrete")))


# For backwards compatibility with iris
ALL_COMPLEX_WAV = available_dt_filters()
ALL_FIRST_STAGE = available_first_stage_filters()


def dtcwt(data, first_stage, wavelet, mode="constant", level=None, axis=-1):
    """
    1D Dual-tree complex wavelet transform [1]_ along an axis.

    Parameters
    ----------
    data: array_like
        Input data. Can be of any shape, but the transform can only be applied in 1D (i.e. along a single axis).
        The length along ``axis`` must be even.
    first_stage : str
        Wavelet to use for the first stage. See :func:`skued.available_first_stage_filters` for a list of suitable arguments
    wavelet : str
        Wavelet to use in stages > 1. Must be appropriate for the dual-tree complex wavelet transform.
        See :func:`skued.available_dt_filters` for possible values.
    mode : str, optional
        Signal extension mode, see :data:`pywt.Modes`.
    level : int or None, optional
        Decomposition level (must be >= 0). If level is None (default) then it
        will be calculated using the ``dt_max_level`` function.
    axis : int, optional
        Axis over which to compute the transform. Default is -1.

    Returns
    -------
    [cA_n, cD_n, cD_n-1, ..., cD2, cD1] : list
        Ordered list of coefficients arrays
        where `n` denotes the level of decomposition. The first element
        (`cA_n`) of the result is approximation coefficients array and the
        following elements (`cD_n` - `cD_1`) are details coefficients arrays.
        Arrays have the same number of dimensions as the input.

    Raises
    ------
    ValueError: Raised if axis argument is invalid (e.g. too large) or if the input data
                is not even along the transform direction.

    Notes
    -----
    The implementation uses two tricks [1]_:

        * Different first-stage wavelet: The first level of the dual-tree
          complex wavelet transform involves a combo of shifted wavelets.
        * Swapping of filters at each stage: At each level > 1, the filters
          (separated into real and imaginary trees) are swapped.

    References
    ----------
    .. [1] Selesnick, I. W. et al. 'The Dual-tree Complex Wavelet Transform', IEEE Signal Processing
           Magazine pp. 123 - 151, November 2005.
    """
    data = np.asarray(data, dtype=float) / np.sqrt(2)

    if level is None:
        level = dt_max_level(
            data=data, first_stage=first_stage, wavelet=wavelet, axis=axis
        )
    elif level == 0:
        return [data]

    # Check axis bounds
    if axis > data.ndim - 1:
        raise ValueError(
            f"Input array has {data.ndim} dimensions, but input axis is {axis}"
        )
    elif data.shape[axis] % 2:
        raise ValueError(
            f"Input array has shape {data.shape[axis]} along transform direction \
                          (axis = {axis}). Even length is required."
        )

    real_wavelet, imag_wavelet = dualtree_wavelet(wavelet)
    real_first, imag_first = dt_first_stage(first_stage)

    real_coeffs = _single_tree_analysis_1d(
        data=data,
        first_stage=real_first,
        wavelet=(real_wavelet, imag_wavelet),
        level=level,
        mode=mode,
        axis=axis,
    )
    imag_coeffs = _single_tree_analysis_1d(
        data=data,
        first_stage=imag_first,
        wavelet=(imag_wavelet, real_wavelet),
        level=level,
        mode=mode,
        axis=axis,
    )

    # Due to the unpredictability of wavelet coefficient sizes,
    # the following operation cannot be done using ndarrays
    return [real + 1j * imag for real, imag in zip(real_coeffs, imag_coeffs)]


def idtcwt(coeffs, first_stage, wavelet, mode="constant", axis=-1):
    """
    1D Inverse dual-tree complex wavelet transform [1]_ along an axis.

    Parameters
    ----------
    coeffs : array_like
        Coefficients list [cAn, cDn, cDn-1, ..., cD2, cD1]
    first_stage : str
        Wavelet to use for the first stage. See :func:`skued.available_first_stage_filters` for a list of suitable arguments
    wavelet : str
        Wavelet to use in stages > 1. Must be appropriate for the dual-tree complex wavelet transform.
        See :func:`skued.available_dt_filters` for possible values.
    mode : str, optional
        Signal extension mode, see :data:`pywt.Modes`.
    axis : int, optional
        Axis over which to compute the inverse transform.

    Returns
    -------
    reconstructed : ndarray
        Reconstructed array

    Raises
    ------
    ValueError: If the input coefficients are too few

    References
    ----------
    .. [1] Selesnick, I. W. et al. 'The Dual-tree Complex Wavelet Transform', IEEE Signal Processing
           Magazine pp. 123 - 151, November 2005.
    """
    if len(coeffs) < 1:
        raise ValueError(
            f"Coefficient list too short with {len(coeffs)} elements (minimum 1 array required)."
        )
    elif len(coeffs) == 1:  # level 0 inverse transform
        return np.sqrt(2) * coeffs[0]
    else:
        real_wavelet, imag_wavelet = dualtree_wavelet(wavelet)
        real_first, imag_first = dt_first_stage(first_stage)

        # reconstruction of real part is first, then
        # reconstruction of imaginary part
        # Large performance boost by adding imaginary part in-place
        inverse = _single_tree_synthesis_1d(
            coeffs=[coeff.real for coeff in coeffs],
            first_stage=real_first,
            wavelet=(real_wavelet, imag_wavelet),
            mode=mode,
            axis=axis,
        )
        inverse += _single_tree_synthesis_1d(
            coeffs=[coeff.imag for coeff in coeffs],
            first_stage=imag_first,
            wavelet=(imag_wavelet, real_wavelet),
            mode=mode,
            axis=axis,
        )

        return np.sqrt(2) * inverse / 2


def dt_max_level(data, first_stage, wavelet, axis=-1):
    """
    Returns the maximum decomposition level possible from the dual-tree complex wavelet transform.

    Parameters
    ----------
    data : ndarray
        Input data. Can be of any dimension.
    first_stage : str
        Wavelet to use for the first stage. See :func:`skued.available_first_stage_filters` for a list of suitable arguments
    wavelet : str
        Wavelet to use in stages > 1. Must be appropriate for the dual-tree complex wavelet transform.
        See :func:`skued.available_dt_filters` for possible values.
    axis : int, optional
        Axis over which to compute the transform. Default is -1

    Returns
    -------
    max_level : int
        Maximum decomposition level.
    """
    real_wavelet, imag_wavelet = dualtree_wavelet(wavelet)
    return dwt_max_level(
        data_len=data.shape[axis],
        filter_len=max([real_wavelet.dec_len, imag_wavelet.dec_len]),
    )


def _normalize_size_axis(approx, detail, axis):
    """
    Adjust the approximate coefficients' array size to that of the detail coefficients' array.

    Parameters
    ----------
    approx : ndarray
    detail: ndarray
    axis : int

    Returns
    -------
    ndarray
        Same shape as detail input.
    """
    if approx.shape[axis] == detail.shape[axis]:
        return approx
    # Swap axes to bring the specific axis to front, truncate array, and re-swap.
    # This is an extension of the 1D case:
    # >>> approx = approx[:-1]
    return np.swapaxes(np.swapaxes(approx, axis, 0)[:-1], 0, axis)


def _single_tree_analysis_1d(data, first_stage, wavelet, level, mode, axis):
    """
    Single tree of the forward dual-tree complex wavelet transform.

    Parameters
    ----------
    data : ndarray, ndim 1
    first_stage : Wavelet object
    wavelet : 2-tuple of Wavelet object
    level : int
    mode : str
    axis : int

    Returns
    -------
    [cA_n, cD_n, cD_n-1, ..., cD2, cD1] : list
        Ordered list of coefficients arrays
        where `n` denotes the level of decomposition. The first element
        (`cA_n`) of the result is approximation coefficients array and the
        following elements (`cD_n` - `cD_1`) are details coefficients arrays.
    """
    approx, first_detail = dwt(data=data, wavelet=first_stage, mode=mode, axis=axis)
    # Use of a deque vs. a list is because deque.appendleft is O(1)
    coeffs_list = deque([first_detail])
    for i, wav in zip(range(level - 1), cycle(wavelet)):
        approx, detail = dwt(data=approx, wavelet=wav, mode=mode, axis=axis)
        coeffs_list.appendleft(detail)

    # Format list ot be compatible to PyWavelet's format. See pywt.wavedec source.
    coeffs_list.appendleft(approx)
    return coeffs_list


def _single_tree_synthesis_1d(coeffs, first_stage, wavelet, mode, axis):
    """
    Single tree of the inverse dual-tree complex wavelet transform.

    Parameters
    ----------
    coeffs : list
    first_stage : Wavelet object
    wavelet : 2-tuple of Wavelet objects
    mode : str
    axis : int

    Returns
    -------
    reconstructed : ndarray, ndim 1
    """
    # Determine the level except first stage:
    # The order of wavelets depends on whether
    # the level is even or odd.
    level = len(coeffs) - 1
    approx, detail_coeffs, first_stage_detail = coeffs[0], coeffs[1:-1], coeffs[-1]

    # Set the right alternating order between real and imag wavelet
    if level % 2 == 1:
        wavelet = reversed(wavelet)

    # In the case of level = 1, coeffs[1:-1] is an empty list. Then, the following
    # loop is not run since zip() iterates as long as the shortest iterable.
    for detail, wav in zip(detail_coeffs, cycle(wavelet)):
        approx = _normalize_size_axis(approx, detail, axis=axis)
        approx = idwt(cA=approx, cD=detail, wavelet=wav, mode=mode, axis=axis)

    approx = _normalize_size_axis(approx, first_stage_detail, axis=axis)
    return idwt(
        cA=approx, cD=first_stage_detail, wavelet=first_stage, mode=mode, axis=axis
    )


@lru_cache(maxsize=len(available_dt_filters()))
def dualtree_wavelet(name):
    """
    Returns a complex wavelet suitable for dual-tree cwt from a name.

    Parameters
    ----------
    name : str

    Returns
    -------
    real, imag : pywt.Wavelet objects.

    Raises
    ------
    ValueError : If illegal wavelet name.
    """
    filters = ("h0a", "h0b", "g0a", "g0b", "h1a", "h1b", "g1a", "g1b")

    filename = DATADIR / (name + ".npz")
    with np.load(filename) as mat:
        try:
            (
                dec_real_low,
                dec_imag_low,
                rec_real_low,
                rec_imag_low,
                dec_real_high,
                dec_imag_high,
                rec_real_high,
                rec_imag_high,
            ) = tuple([mat[k].flatten() for k in filters])
        except KeyError:
            raise ValueError(f"Wavelet does not define ({filters}) coefficients")

    real_filter_bank = [dec_real_low, dec_real_high, rec_real_low, rec_real_high]
    imag_filter_bank = [dec_imag_low, dec_imag_high, rec_imag_low, rec_imag_high]

    return (
        Wavelet(name="real:" + name, filter_bank=real_filter_bank),
        Wavelet(name="imag:" + name, filter_bank=imag_filter_bank),
    )


@lru_cache(maxsize=len(available_first_stage_filters()))
def dt_first_stage(wavelet):
    """
    Returns two wavelets to be used in the dual-tree complex wavelet transform, at the first stage.

    Parameters
    ----------
    wavelet : str or Wavelet

    Return
    ------
    wav1, wav2 : Wavelet objects

    Raises
    ------
    ValueError : If invalid first stage wavelet.
    """
    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)

    if wavelet.name not in available_first_stage_filters():
        raise ValueError(f"{wavelet.name} is an invalid first stage wavelet.")

    # extend filter bank with zeros
    filter_bank = [np.array(f, copy=True) for f in wavelet.filter_bank]
    for i, filt in enumerate(filter_bank):
        extended = np.zeros(shape=(filt.shape[0] + 2,), dtype=float)
        extended[1:-1] = filt
        filter_bank[i] = extended

    # Shift deconstruction filters to one side, and reconstruction
    # to the other side
    shifted_fb = [np.array(f, copy=True) for f in filter_bank]
    for i, filt in enumerate(shifted_fb[:2]):  # Deconstruction filters
        shifted_fb[i] = np.roll(filt, 1)
    for i, filt in enumerate(shifted_fb[2:], start=2):  # Reconstruction filters
        shifted_fb[i] = np.roll(filt, -1)

    return (
        Wavelet(name=wavelet.name, filter_bank=filter_bank),
        Wavelet(name=wavelet.name, filter_bank=shifted_fb),
    )
